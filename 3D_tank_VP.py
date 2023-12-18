import pdb
import time
import numpy as np
import os.path
from firedrake import *
from savings import save_README
from firedrake.petsc import PETSc
from FIAT.quadrature import GaussLobattoLegendreQuadratureLineRule
from FIAT.reference_element import UFCInterval

# Input the test case you are going to run below:
case = 'TC4' # TC1/TC3/TC4/TCU
assert case.upper() in ['TC1', 'TC3', 'TC4', 'TCU'], "Incorrect input!"

start_time = time.perf_counter()

"""
    ****************************************
    *            Load Settings             *
    **************************************** """
PETSc.Sys.Print('Setting up test case %3s...' % case)

if case=='TC1':
    from settings_TC1 import *
    input_data, scheme, dim, save_path, bottom, FWF, save_pvd = test_case()
    k, H0, xb, sb, H_expr, Hend, Lx, Ly, Lw, res_x, res_y, n_z = domain(bottom)
    g, w, Tw, gamma, t_stop, Aa, Bb, WM_expr, dWM_dt_expr, dWM_dy_expr, h_ex_expr, phis_ex_expr, phii_ex_expr = wavemaker(k, dim, H0, Ly, Lw, input_data)
    T0, t, dt, Tend, dt_save = set_time(Lx, res_x, Ly, res_y, H0, n_z, Tw)
else:
    if case=='TC2': 
        from settings_TC2 import *
    elif case=='TC3': 
        from settings_TC3 import *
    elif case=='TC4':
        from settings_TC4 import *
    else:
        from settings_User import *
    input_data, scheme, dim, save_path, bottom, FWF, save_pvd = test_case()
    H0, xb, sb, H_expr, Hend, Lx, Ly, Lw, res_x, res_y, n_z = domain(bottom)
    g, lamb, k, w, Tw, gamma, t_stop, WM_expr, dWM_dt_expr, dWM_dy_expr = wavemaker(dim, H0, Ly, Lw, input_data)
    T0, t, dt, Tend, dt_save = set_time(Tw)

one_ver_ele = True # Whether there's only one element in the vertical direction of the 3D extruded mesh
hatphi_one = False # Whether hat_phi(z) is set to be a constant one.

PETSc.Sys.Print('...settings loaded!')

# Create the directory to store output files
try:
    os.makedirs(save_path, exist_ok = True)
    PETSc.Sys.Print("Directory created successfully!")
except OSError as error:
    PETSc.Sys.Print("Directory can not be created!")

"""
    ****************************************************
    *               Definition of the mesh             *
    **************************************************** """

PETSc.Sys.Print('Creation of the mesh across %d processes...' % COMM_WORLD.size)

#_________________ Vertical discretization ________________#
Nz = n_z+1         # Number of layers in one vertical element

#________________ Horizontal discretization _______________#
Nx = round(Lx/res_x)    # Number of elements in x (round to the nearest integer)
Ny = round(Ly/res_y)    # Number of elements in y

#___________________________ Mesh _________________________#
if dim=='3D':
    hor_mesh = RectangleMesh(Nx, Ny, Lx, Ly, quadrilateral=True)
elif dim=='2D':
    hor_mesh = IntervalMesh(Nx, Lx)
    #The left hand boundary point has boundary marker 1, while the right hand point has marker 2.
# see https://www.firedrakeproject.org/extruded-meshes.html
if one_ver_ele:
    mesh = ExtrudedMesh(hor_mesh, 1, layer_height=H0, extrusion_type='uniform')
    PETSc.Sys.Print('-> Only 1 element in the vertical direction.')
else:
    mesh = ExtrudedMesh(hor_mesh, n_z, layer_height=H0/n_z, extrusion_type='uniform')
    PETSc.Sys.Print(f'-> {n_z} elements in the vertical direction.')

x = SpatialCoordinate(mesh)

PETSc.Sys.Print('...extruded mesh created!')

"""
    *************************************************
    *          Definition of the functions          *
    ************************************************* """
PETSc.Sys.Print('Definition of the function...')

#__________________Define function space___________________#
nCG = 1  # degree of the function space in the horizontal direction
if one_ver_ele: # determine the degree of the function space in the vertical direction
    nCG_v = n_z
else:
    nCG_v = 1
    
V   = FunctionSpace(mesh, 'CG', nCG, vfamily='CG', vdegree=nCG_v)
V_h = FunctionSpace(mesh, 'CG', nCG, vfamily='R', vdegree=0)
V_mixed = V_h * V_h * V

#_________________________Unknows__________________________#
# At time step t^n
h_n0   = Function(V_h)      
psi_n0 = Function(V_h)                           
phi_n0 = Function(V) # full velocity potential, only for the purpose of output                    

# At time step t^{n+1}
h_n1   = Function(V_h)                      
psi_n1 = Function(V_h)

# In the middle t^{n+1/2}
solutions = Function(V_mixed)
psi_mp, h_mp, varphi_mp = split(solutions)

#___________Parition of the velocity potentials____________#
hat_phi = Function(V)
hat_phi_z = Function(V)

#____________________Wave tank geometry____________________# 
b = Function(V_h) # bathymetry b(x,y)
H = Function(V_h) # rest water depth H(x,y)

#________________________Wavemaker_________________________#
WM_n0   = Function(V_h)  
WM_n1   = Function(V_h)
WM_half = Function(V_h)
dWM_half_dy = Function(V_h)
dWM_half_dt = Function(V_h)

#___________________Auxiliary functions____________________#
X = Function(V_h)
W_n0   = Function(V_h)
W_n1   = Function(V_h)
W_half = Function(V_h)
U_half = Function(V_h)

#______________________Test functions______________________#
w_t = TestFunction(V_mixed)   
delta_psi, delta_h, delta_varphi = split(w_t)

PETSc.Sys.Print('...functions created!')

"""
    ***********************************************************************************
    *                          Initialisation of the Functions                        *
    ***********************************************************************************"""
PETSc.Sys.Print('Initalisation of the functions...')

#____________________Wave tank geometry____________________#
H_expr(H,x)     # Water depth at rest H(x,y)
b.assign(H0-H)  # Bathymetry b(x,y)

#________________________Wavemaker_________________________#
t_half = t + 0.5*dt
t_n1 = t + dt

if input_data=="measurements":
    with open('202002/PistonMotion.dat','r') as data_f1:
        t_wm, wm_data = np.loadtxt(data_f1, usecols=(0,1), unpack=True) # time and wavemaker position
    with open('202002/PistonVelocity.dat','r') as data_f2:
        t_vel, vel_data = np.loadtxt(data_f2, usecols= (0,1), unpack=True) # time and wavemaker velocity
    WM_expr(WM_n0, x, t, t_wm, wm_data)
    WM_expr(WM_n1, x, t_n1, t_wm, wm_data)
    WM_expr(WM_half, x, t_half, t_wm, wm_data)
    dWM_dt_expr(dWM_half_dt, x, t_half, t_vel, vel_data)
    dWM_dy_expr(dWM_half_dy)
else:
    WM_expr(WM_n0,  x, t, t_stop)
    WM_expr(WM_n1,  x, t_n1, t_stop)
    WM_expr(WM_half,x, t_half, t_stop)               
    dWM_dt_expr(dWM_half_dt, x, t_half, t_stop) 
    dWM_dy_expr(dWM_half_dy, x, t_half, t_stop)        

#___________________Auxiliary functions____________________#
X.interpolate(x[0]-Lw)
U_half.interpolate(X*dWM_half_dy)
W_n0.assign(Lw-WM_n0)
W_n1.assign(Lw-WM_n1)
W_half.assign(Lw-WM_half)

#___________Parition of the velocity potentials____________#
if hatphi_one:
    hat_phi   = Constant(1.0)
    hat_phi_z = Constant(0.0)
    PETSc.Sys.Print('-> hatphi = 1')
    # hat_phi.interpolate(x[1]/H0)
    # hat_phi_z = Constant(1/H0)
    # PETSc.Sys.Print('-> hatphi = z/H0') # can improve the accuracy at the bottom (z=0)
else: # Construct hat_phi(z) as Lagrange polynomial based on GLL points
    PETSc.Sys.Print(f'-> hatphi = GLL of order {n_z}')
    fiat_rule = GaussLobattoLegendreQuadratureLineRule(UFCInterval(), Nz)
    z_k = H0*fiat_rule.get_points() # np.array, len(z_k)=Nz
    if dim =='3D':
        hat_phi_expr = product( (x[2]-z_k.item(i))/(H0-z_k.item(i)) for i in range(n_z) )
        hat_phi.interpolate(hat_phi_expr)
        hat_phi_z.interpolate(hat_phi.dx(2))
    elif dim=='2D':
        hat_phi_expr = product( (x[1]-z_k.item(i))/(H0-z_k.item(i)) for i in range(n_z) )
        hat_phi.interpolate(hat_phi_expr)
        hat_phi_z.interpolate(hat_phi.dx(1))

#__________________Seeding the solutions___________________#
if case=='TC1':
    h_ex_expr(h_n0,x,t)
    phis_ex_expr(psi_n0,x,t)
    # For TC1, it's necessary to set initial values for psi_mp and h_mp
    solutions.sub(0).assign(psi_n0)
    solutions.sub(1).assign(H)
    # Will the results change if varphi_mp is also initialised? No!
    #phii_ex_expr(phi_n0,x,t)
    #solutions.sub(2).interpolate(phi_n0-psi_n0*hat_phi)
else:
    h_n0.assign(H)  # h(x,y;t=0) = H(x)
    solutions.sub(1).assign(H)

PETSc.Sys.Print('...functions initialised!')

"""
    **************************************************************************************
    *                                 Define the solvers                                 *
    ************************************************************************************** """
PETSc.Sys.Print('Initialisation of the solvers...')

''' # 3D VP with flat bottom and no wavemaker motion
    VP_mp= (-psi_mp*(h_n1-h_n0)/dt + h_mp*(psi_n1-psi_n0)/dt )*ds_t\
          +( (h_mp/(2*H0)) * ( hat_phi*psi_mp.dx(0)+varphi_mp.dx(0)- (x[2]*h_mp.dx(0)/h_mp)*(psi_mp*hat_phi_z+varphi_mp.dx(2)) )**2 \
            +(h_mp/(2*H0)) * ( hat_phi*psi_mp.dx(1)+varphi_mp.dx(1)- (x[2]*h_mp.dx(1)/h_mp)*(psi_mp*hat_phi_z+varphi_mp.dx(2)) )**2 \
            +(H0/(2*h_mp)) * ( psi_mp*hat_phi_z + varphi_mp.dx(2) )**2 \
            + (g*h_mp/H0)  * ( x[2]*h_mp/H0 - H0 )  )*dx
'''

if dim=='2D':
    VP_mp= (H0*( psi_mp*W_half*(h_n1-h_n0)/dt - h_mp*(psi_n1*W_n1-psi_n0*W_n0)/dt + psi_mp*X*dWM_half_dt*h_mp.dx(0) ))*ds_t \
        -( (Lw*Lw*h_mp/(2*W_half)) * ( hat_phi*psi_mp.dx(0)+varphi_mp.dx(0) - ((H0*b.dx(0)+x[1]*h_mp.dx(0))/h_mp)*(psi_mp*hat_phi_z+varphi_mp.dx(1)) )**2 \
         + (H0*H0*W_half/(2*h_mp)) * ( psi_mp*hat_phi_z+varphi_mp.dx(1) )**2 \
          + W_half*g*h_mp*( (x[1]*h_mp/H0) - H ) )*dx \
        - ( Lw*dWM_half_dt*h_mp*(psi_mp*hat_phi+varphi_mp) )*ds_v(1)
elif dim=='3D':
    VP_mp= (H0*( psi_mp*W_half*(h_n1-h_n0)/dt - h_mp*(psi_n1*W_n1-psi_n0*W_n0)/dt + psi_mp*X*dWM_half_dt*h_mp.dx(0) ))*ds_t \
        -( (Lw*Lw*h_mp/(2*W_half)) * ( hat_phi*psi_mp.dx(0)+varphi_mp.dx(0) - ((H0*b.dx(0)+x[2]*h_mp.dx(0))/h_mp)*(psi_mp*hat_phi_z+varphi_mp.dx(2)) )**2 \
               + (0.5*W_half*h_mp) * ( hat_phi*psi_mp.dx(1)+varphi_mp.dx(1) - ((H0*b.dx(1)+x[2]*h_mp.dx(1))/h_mp)*(psi_mp*hat_phi_z+varphi_mp.dx(2)) \
                     +(U_half/W_half)*(hat_phi*psi_mp.dx(0)+varphi_mp.dx(0) - ((H0*b.dx(0)+x[2]*h_mp.dx(0))/h_mp)*(psi_mp*hat_phi_z+varphi_mp.dx(2))) )**2 \
         + (H0*H0*W_half/(2*h_mp)) * ( psi_mp*hat_phi_z+varphi_mp.dx(2) )**2 \
          + W_half*g*h_mp*( (x[2]*h_mp/H0) - H ) )*dx \
        - ( Lw*dWM_half_dt*h_mp*(psi_mp*hat_phi+varphi_mp) )*ds_v(1)

#ds_t is used to denote an integral over the top surface of the mesh
#ds_v is used to denote an integral over side facets of the mesh.

psi_expr = derivative(VP_mp, psi_mp, du=delta_psi) 
psi_expr = replace(psi_expr, {psi_n1: 2.0*psi_mp-psi_n0})
psi_expr = replace(psi_expr, {h_n1: 2.0*h_mp-h_n0})
# https://www.firedrakeproject.org/firedrake.html#firedrake.ufl_expr.derivative

h_expr = derivative(VP_mp, h_mp, du=delta_h)
h_expr = replace(h_expr, {psi_n1: 2.0*psi_mp-psi_n0})
h_expr = replace(h_expr, {h_n1: 2.0*h_mp-h_n0})

phi_expr = derivative(VP_mp, varphi_mp, du=delta_varphi)
phi_expr = replace(phi_expr, {psi_n1: 2.0*psi_mp-psi_n0})
phi_expr = replace(phi_expr, {h_n1: 2.0*h_mp-h_n0})

# BC
BC_varphi = DirichletBC(V_mixed.sub(2), Constant(0.0), 'top')
# top, to set a boundary condition on the top surface.
'''
sol_para = {'ksp_type': 'gmres', 
            'pc_type': 'python', 
            'pc_python_type': 'firedrake.ASMStarPC', 
            'star_construct_dim': 2,
            'star_sub_sub_pc_type': 'lu', 
            'sub_sub_pc_factor_mat_ordering_type': 'rcm'}
'''
sol_para = {'ksp_converged_reason': None}

WF_expr =  psi_expr + h_expr + phi_expr 
problem_mp = NonlinearVariationalProblem(WF_expr, solutions, bcs = BC_varphi)
solver_mp = NonlinearVariationalSolver(problem_mp, solver_parameters=sol_para)

PETSc.Sys.Print('...solver initialised!')

"""
    *************************************************************
    *                        Saving Files                       *
    ************************************************************* """
readme_file = os.path.join(save_path, 'readme.txt')
energy_file = os.path.join(save_path, 'energy.csv')
energy_data = np.empty((0,4)) # create an empty 2d array for saving E(t)

if case=='TC4':
    probe_file = os.path.join(save_path, 'probes.csv')
    probe_data = np.empty((0,9))

""" *********************************************************************************
    *                                   Time loop                                   *
    ********************************************************************************* """
PETSc.Sys.Print('Update of the solutions:')

t_save = t+dt_save # t_save = t + dt*m, output starts at the m-th time step
smallfac = 10.0**(-10)

before_it = time.perf_counter()-start_time # running time from start until this line
#pdb.set_trace()

while t<=Tend+smallfac: 
    
    #_____________________Call the solver______________________#
    solver_mp.solve()
    psi_mp, h_mp, varphi_mp = solutions.subfunctions

    #______________________Next time step______________________#
    t = t + dt
    psi_n0.assign(2*psi_mp-psi_n0)
    h_n0.assign(2*h_mp-h_n0)
    if case=='TC1':
        phi_n0.interpolate(psi_n0*hat_phi+varphi_mp) # approximation used! Only for the purpose of output

    #_____________________Save the results_____________________#
    if t_save-smallfac < t: 
        
        progress = format(100*t/Tend, '.3f')+' %'
        tt = format(t, '.4f')
        PETSc.Sys.Print('t= %s, Progress: %s' % (tt, progress))

        # Calculate energy using the time-discretised VP
        if dim=='2D':
            kin_energy=assemble( ( (Lw*Lw*h_mp/(2*W_half)) * ( hat_phi*psi_mp.dx(0)+varphi_mp.dx(0) - ((H0*b.dx(0)+x[1]*h_mp.dx(0))/h_mp)*(psi_mp*hat_phi_z+varphi_mp.dx(1)) )**2 \
                                  +(H0*H0*W_half/(2*h_mp)) * ( psi_mp*hat_phi_z+varphi_mp.dx(1) )**2 )*dx )
            pot_energy=assemble( (W_half*g*h_mp*( (x[1]*h_mp/H0) - H ))*dx )
        elif dim=='3D':
            kin_energy=assemble( ((Lw*Lw*h_mp/(2*W_half)) * ( hat_phi*psi_mp.dx(0)+varphi_mp.dx(0) - ((H0*b.dx(0)+x[2]*h_mp.dx(0))/h_mp)*(psi_mp*hat_phi_z+varphi_mp.dx(2)) )**2 \
                                      + (0.5*W_half*h_mp) * ( hat_phi*psi_mp.dx(1)+varphi_mp.dx(1) - ((H0*b.dx(1)+x[2]*h_mp.dx(1))/h_mp)*(psi_mp*hat_phi_z+varphi_mp.dx(2)) \
                                         +(U_half/W_half) * ( hat_phi*psi_mp.dx(0)+varphi_mp.dx(0) - ((H0*b.dx(0)+x[2]*h_mp.dx(0))/h_mp)*(psi_mp*hat_phi_z+varphi_mp.dx(2))) )**2 \
                                + (H0*H0*W_half/(2*h_mp)) * ( psi_mp*hat_phi_z + varphi_mp.dx(2) )**2 )*dx )
            pot_energy=assemble( (W_half*g*h_mp*( (x[2]*h_mp/H0) - H ))*dx )
        
        tot_energy = kin_energy + pot_energy

        energy_data = np.append(energy_data,np.array([[t,kin_energy,pot_energy,tot_energy]]),axis=0)

        if case=='TC1':
            # Store the field data to a numpy array
            field_data = np.empty((Nx+1,5))
            for row_i in range(Nx+1):
                x_i = row_i*(Lx/Nx) # x coordinates
                field_data[row_i,:]=np.array([x_i, h_n0.at([x_i,0]), psi_n0.at([x_i,0]), phi_n0.at([x_i,H0]), phi_n0.at([x_i,0])])

            # Save the array to a binary file (.npy)
            op_file = os.path.join(save_path, tt+'.npy')
            with open(op_file,'wb') as f:
                np.save(f, field_data)

        if case=='TC4':
            probe_data = np.append(probe_data,\
                np.array([[t,*h_n0.at([10,0],[20,0],[40,0],[49.5,0],[50,0],[54,0]),WM_n1.at([0,0]),dWM_half_dt.at([0,0])]]),axis=0)

        t_save+=dt_save

    
    #___________________Update the wavemaker___________________#
    t_half = t + 0.5*dt
    t_n1   = t + dt
    
    if input_data=="created":
        WM_expr(WM_n0,   x, t,      t_stop)
        WM_expr(WM_half, x, t_half, t_stop)
        WM_expr(WM_n1,   x, t_n1,   t_stop)
        dWM_dt_expr(dWM_half_dt, x, t_half, t_stop) 
        dWM_dy_expr(dWM_half_dy, x, t_half, t_stop)
    else: # input_data=='measurements'
        WM_expr(WM_n0, x, t, t_wm, wm_data)
        WM_expr(WM_n1, x, t_n1, t_wm, wm_data)
        WM_expr(WM_half, x, t_half, t_wm, wm_data)
        dWM_dt_expr(dWM_half_dt, x, t_half, t_vel, vel_data)

    W_n0.assign(Lw-WM_n0)
    W_n1.assign(Lw-WM_n1)
    W_half.assign(Lw-WM_half)
    U_half.interpolate(X*dWM_half_dy)

with open(energy_file,'w') as e_f:
    np.savetxt(e_f, energy_data, fmt='%10.4f '+' %.18e '*3)

if case=='TC4':
    with open(probe_file,'w') as p_f:
        np.savetxt(p_f, probe_data, fmt='%10.4f'*9)

comp_time = time.perf_counter()-start_time
jours = int(comp_time/(24*3600))
heures = int((comp_time-jours*24*3600)/3600)
minutes = int((comp_time-jours*24*3600-heures*3600)/60)
secondes = comp_time -jours*24*3600-heures*3600 - minutes*60
with open(readme_file,'w') as info:
    save_README(info, dim, Lx, Ly, H0, xb, sb, Nx, Ny, n_z, gamma, Tw, w, t_stop, Lw, scheme, dt, t,\
                jours, heures, minutes, secondes, comp_time, before_it)
