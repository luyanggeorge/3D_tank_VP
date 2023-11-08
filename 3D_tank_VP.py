import pdb
import time
import numpy as np
import os.path
from scipy.interpolate import interp1d
from vertical_discr_full import *
from savings import *
from firedrake import *
from firedrake.petsc import PETSc

# Input the test case you are going to run below:
case = 'TC3' # TC3/TC4/TCU
assert case.upper() in ['TC3', 'TC4', 'TCU'], "Incorrect input!"

start_time = time.perf_counter()

"""
    ****************************************
    *            Load Settings             *
    **************************************** """
PETSc.Sys.Print('Setting up test case %3s...' % case)

if case=='TC3': 
    from settings_TC3 import *
elif case=='TC4':
    from settings_TC4 import *
else:
    from settings_User import *

input_data, scheme, dim, save_path, bottom, FWF, save_pvd = test_case()
H0, xb, sb, H_expr, Hend, Lx, Ly, Lw, res_x, res_y, n_z = domain(bottom)
g, lamb, k, w, Tw, gamma, t_stop, WM_expr, dWM_dt_expr, dWM_dy_expr = wavemaker(dim, H0, Ly, Lw, input_data)
T0, t, dt, Tend, dt_save = set_time(Tw)

PETSc.Sys.Print('...settings loaded!')

    
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
# see https://www.firedrakeproject.org/extruded-meshes.html
hor_mesh = RectangleMesh(Nx, Ny, Lx, Ly, quadrilateral=True)
mesh = ExtrudedMesh(hor_mesh, n_z, layer_height=H0/n_z, extrusion_type='uniform', name='mesh')
x = SpatialCoordinate(mesh)

PETSc.Sys.Print('...mesh created!')

"""
    *************************************************
    *          Definition of the functions          *
    ************************************************* """
PETSc.Sys.Print('Definition of the function...')

#__________________Define function space___________________#
nCG = 1  # degree of the function space in the horizontal direction
nCG_v= 1 # degree of the function space in the vertical direction
V   = FunctionSpace(mesh, 'CG', nCG, vfamily='CG', vdegree=nCG_v)
V_h = FunctionSpace(mesh, 'CG', nCG, vfamily='R', vdegree=0)
V_mixed = V_h * V_h * V

#_________________________Unknows__________________________#
# At time step t^n
h_n0   = Function(V_h)      
psi_n0 = Function(V_h)                           
#varphi_n0 = Function(V)                     

# At time step t^{n+1}
h_n1   = Function(V_h)                      
psi_n1 = Function(V_h)    
#varphi_n1 = Function(V)

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
H.assign(H_expr(V_h,x))    # Depth at rest H(x)
b.assign(H0-H)             # Beach b(x,y)

#________________________Wavemaker_________________________#
t_half = t + 0.5*dt
t_n1 = t + dt

if input_data=="measurements":
    # to be changed for TC4
    '''
    t_data, wm_data, wm_vel_data = load_wavemaker(dt)
    wm_inter     = interp1d(t_data,wm_data)
    wm_vel_inter = interp1d(t_data,wm_vel_data)
    
    WM_expr(WM, x, t, wm_inter)
    dWM_dt_expr(dWM_dt, x, t, wm_vel_inter)
    dWM_dy_expr(dWM_dy)
    '''
else:
    WM_expr(WM_n0,  x, t)
    WM_expr(WM_n1,  x, t_n1)
    WM_expr(WM_half,x, t_half)               
    dWM_dt_expr(dWM_half_dt,x,t_half) 
    dWM_dy_expr(dWM_half_dy,x,t_half)        

#___________________Auxiliary functions____________________#
X.interpolate(x[0]-Lw)
U_half.interpolate(X*dWM_half_dy)
W_n0.assign(Lw-WM_n0)
W_n1.assign(Lw-WM_n1)
W_half.assign(Lw-WM_half)


#__________________Seeding the solutions___________________#
h_n0.assign(H)  # h(x,y;t=0) = H(x)

#___________Parition of the velocity potentials____________#
hat_phi.assign(Constant(1.0))
hat_phi_z.assign(Constant(0.0))

PETSc.Sys.Print('...functions initialised!')

"""
    **************************************************************************************
    *                                 Define the solvers                                 *
    ************************************************************************************** """
PETSc.Sys.Print('Initialisation of the solvers...')

VP_mp =  (H0*( psi_mp*W_half*(h_n1-h_n0)/dt - h_mp*(psi_n1*W_n1-psi_n0*W_n0)/dt + psi_mp*X*dWM_half_dt*h_mp.dx(0) ))*ds_t \
        -( (Lw*Lw*h_mp/(2*W_half)) * ( hat_phi*psi_mp.dx(0)+varphi_mp.dx(0) - ((H0*b.dx(0)+x[2]*h_mp.dx(0))/h_mp)*(psi_mp*hat_phi_z+varphi_mp.dx(2)) )**2 \
               + (0.5*W_half*h_mp) * ( hat_phi*psi_mp.dx(1)+varphi_mp.dx(1) - ((H0*b.dx(1)+x[2]*h_mp.dx(1))/h_mp)*(psi_mp*hat_phi_z+varphi_mp.dx(2)) \
                     +(U_half/W_half)*(hat_phi*psi_mp.dx(0)+varphi_mp.dx(0) - ((H0*b.dx(0)+x[2]*h_mp.dx(0))/h_mp)*(psi_mp*hat_phi_z+varphi_mp.dx(2))) )**2 \
         + (H0*H0*W_half/(2*h_mp)) * ( psi_mp*hat_phi_z+varphi_mp.dx(2) )**2 \
          + W_half*g*h_mp*( (x[2]*h_mp/H0) - H ) )*dx \
        - ( Lw*dWM_half_dt*h_mp*(psi_mp*hat_phi+varphi_mp) )*ds_v(1)

#ds_t is used to denote an integral over the top surface of the mesh
#ds_v is used to denote an integral over side facets of the mesh.

# Step-1: solve h^(n+1/2) wrt psi^(n+1/2)
psi_expr = derivative(VP_mp, psi_mp, du=delta_psi) 
psi_expr = replace(psi_expr, {psi_n1: 2.0*psi_mp-psi_n0})
psi_expr = replace(psi_expr, {h_n1: 2.0*h_mp-h_n0})

# wave breaking: psi_expr = psi_expr + wave_breaking
# https://www.firedrakeproject.org/firedrake.html#firedrake.ufl_expr.derivative

# Step-2: solve psi^(n+1/2) wrt hmp=h^(n+1/2)
h_expr = derivative(VP_mp, h_mp, du=delta_h)
h_expr = replace(h_expr, {psi_n1: 2.0*psi_mp-psi_n0})
h_expr = replace(h_expr, {h_n1: 2.0*h_mp-h_n0})

# Step-3: wrt varmp=varphi^(n+1/2) solve varmp=varphi^(n+1/2)
phi_expr = derivative(VP_mp, varphi_mp, du=delta_varphi)
phi_expr = replace(phi_expr, {psi_n1: 2.0*psi_mp-psi_n0})
phi_expr = replace(phi_expr, {h_n1: 2.0*h_mp-h_n0})

# BC
BC_varphi = DirichletBC(V_mixed.sub(2), Constant(0), 'top')
# top, to set a boundary condition on the top surface.

sol_para = {'ksp_type': 'gmres', 'pc_type': 'python', 'pc_python_type': 'firedrake.ASMStarPC', 'star_construct_dim': 2,
            'star_sub_sub_pc_type': 'lu', 'sub_sub_pc_factor_mat_ordering_type': 'rcm'}

WF_expr =  psi_expr + h_expr + phi_expr 
#wave breaking: RHS +WF_wave_breaking_psi
problem_mp = NonlinearVariationalProblem(WF_expr, solutions, bcs = BC_varphi)
solver_mp = NonlinearVariationalSolver(problem_mp, solver_parameters=sol_para)

PETSc.Sys.Print('...solver initialised!')

#pdb.set_trace()

"""
    *************************************************************
    *                        Saving Files                       *
    ************************************************************* """





""" *********************************************************************************
    *                                   Time loop                                   *
    ********************************************************************************* """
PETSc.Sys.Print('Update of the solutions:')

t_save = t
before_it = time.perf_counter()-start_time # running time from start until this line
t_aux = t
update_wm = True # It's used to record the time when t>t_stop for the first time.
smallfac = 10.0**(-10)

#pdb.set_trace()

while t<=Tend+smallfac: 
    
    #_____________________Call the solver______________________#
    solver_mp.solve()
    psi_mp, h_mp, varphi_mp = solutions.split()

    #___________________Calculate the energy___________________#
    if t_save-smallfac < t: 
        
        progress = format(100*t/Tend, '.3f')+' %'
        tt = format(t, '.3f')
        PETSc.Sys.Print('t= %s, Progress: %s' % (tt, progress))

        # calculate energy using the time-discretised VP
        kin_energy=assemble( ((Lw*Lw*h_mp/(2*W_half)) * ( hat_phi*psi_mp.dx(0)+varphi_mp.dx(0) - ((H0*b.dx(0)+x[2]*h_mp.dx(0))/h_mp)*(psi_mp*hat_phi_z+varphi_mp.dx(2)) )**2 \
                                  + (0.5*W_half*h_mp) * ( hat_phi*psi_mp.dx(1)+varphi_mp.dx(1) - ((H0*b.dx(1)+x[2]*h_mp.dx(1))/h_mp)*(psi_mp*hat_phi_z+varphi_mp.dx(2)) \
                                     +(U_half/W_half) * ( hat_phi*psi_mp.dx(0)+varphi_mp.dx(0) - ((H0*b.dx(0)+x[2]*h_mp.dx(0))/h_mp)*(psi_mp*hat_phi_z+varphi_mp.dx(2))) )**2 \
                            + (H0*H0*W_half/(2*h_mp)) * ( psi_mp*hat_phi_z + varphi_mp.dx(2) )**2 )*dx)
        pot_energy=assemble( (W_half*g*h_mp*( (x[2]*h_mp/H0) - H ))*dx )
        tot_energy = kin_energy + pot_energy

        # write into a file

        t_save+=dt_save

    #______________________Next time step______________________#
    t = t + dt
    t_half = t + 0.5*dt
    t_n1   = t + dt

    psi_n0.assign(2*psi_mp-psi_n0)
    h_n0.assign(2*h_mp-h_n0)
    
    if input_data=="created":
        if t <= t_stop: # assume t_half <= t_stop                                   

            WM_expr(WM_n0,   x, t)
            WM_expr(WM_half, x, t_half)
            WM_expr(WM_n1,   x, t_n1)
            dWM_dt_expr(dWM_half_dt, x, t_half) 
            dWM_dy_expr(dWM_half_dy, x, t_half)
            
            W_n0.assign(Lw-WM_n0)
            W_n1.assign(Lw-WM_n1)
            W_half.assign(Lw-WM_half)
            U_half.interpolate(X*dWM_half_dy)

            t_aux = t      # yl notes: store the time when the wave maker stops

        elif t>t_stop and update_wm:      # We stop the wavemaker motion;
            update_wm = False
            if scheme=="SV":
                if t_half<=t_stop:
                    t_aux = t_half
            
            WM_expr(WM_n1,x,t_aux)
            dWM_n1_dt.assign(0.0)
            dWM_dy_expr(dWM_n1_dy,x,t_aux)
            Ww_n1.assign(Lw-WM_n1)

            

    #else: # input_data=='measurements'


comp_time = time.perf_counter()-start_time
jours = int(comp_time/(24*3600))
heures = int((comp_time-jours*24*3600)/3600)
minutes = int((comp_time-jours*24*3600-heures*3600)/60)
secondes = comp_time -jours*24*3600-heures*3600 - minutes*60
#save_README(README_file, Lx, Ly, H0, xb, sb, res_x, Nx, Ny, n_z, gamma, Tw, w, t_stop, Lw, scheme, dt, t,\
#            jours, heures, minutes, secondes, comp_time, before_it)
