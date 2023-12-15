# Post-processing code for TC1 using the VP approach
# Validate the numerical results against the standing wave solutions of linear potential-flow eqs
# Calculate the L2 norm for each flame of data

import numpy as np
from numpy.linalg import norm
from matplotlib import pyplot as plt
import os.path

def load_binary(file_path, file_name):
    """load data from binary files (.npy)"""
    file = os.path.join(file_path, file_name)
    with open(file,'rb') as f:
        data_array = np.load(f)
    x_i = data_array[:,0]
    h   = data_array[:,1]
    psi = data_array[:,2]
    phis = data_array[:,3]
    phib = data_array[:,4]
    return x_i, h, psi, phis, phib

#------ User Input ------#
m = 2
Lx = 2*np.pi

Aa = 0.01
Bb = 0.01

data_path ='data_VP/MMP/TC1_test4_2D(one-ele,phiGLL)'
label = r'2D(one-ele, $\hat{\phi}=$GLL)'

interval=1

save_figure=True
comparison =True
data_path2 ='data_VP/MMP/TC1_test1_dt'
label2 = r'3D(mul-ele, $\hat{\phi}=1$)'
#------------------------#

if save_figure:
    figure_name='TC1_L2_comparison(GLL-1).png'
    path=os.path.join(data_path,'figures')
    try: 
        os.makedirs(path, exist_ok = True) 
    except OSError as error: 
        print(error)
    save_path=os.path.join(path, figure_name)

# parameters
k = 2*np.pi*m/Lx
g  = 9.81
kH0 = np.arccosh(0.5*g)
H0 = kH0/k
w = np.sqrt(2*k*np.sinh(k*H0))
amp_eta=np.sqrt(Aa*Aa + Bb*Bb)
amp_phi=amp_eta*(np.exp(k*H0)+np.exp(-k*H0))/w

# get t from the energy file
energy_file = os.path.join(data_path,'energy.csv')
with open(energy_file,'r') as e_f:
    t_steps_full = np.loadtxt(e_f, usecols=0)
t_steps = t_steps_full[::interval]

L2_h, L2_psi, L2_phis, L2_phib = [],[],[],[]

for t in t_steps:
	# numerical results
    tt   = format(t,'.4f')
    fname = tt+'.npy'
    x_n, h_n, psi_n, phis_n, phib_n = load_binary(data_path, fname)

	# exact standing wave solutions
    h_ex = H0 + np.cos(k*x_n) * ( Aa*np.cos(w*t) + Bb*np.sin(w*t) )
    phis_ex = np.cos(k*x_n) * (np.exp(k*H0) + np.exp(-k*H0)) * (-Aa*np.sin(w*t) + Bb*np.cos(w*t))/w
    phib_ex = np.cos(k*x_n) * (np.exp(k*0) + np.exp(-k*0)) * (-Aa*np.sin(w*t) + Bb*np.cos(w*t))/w

    L2_h.append(norm(h_ex-h_n)) 
    L2_psi.append(norm(phis_ex-psi_n))
    L2_phis.append(norm(phis_ex-phis_n))
    L2_phib.append(norm(phib_ex-phib_n))

if comparison:
	# get t from the energy file
    energy_file2 = os.path.join(data_path2,'energy.csv')
    with open(energy_file2,'r') as e_f2:
        t_steps_full2 = np.loadtxt(e_f2, usecols=0)
    t_steps2 = t_steps_full2[::interval]

    L2_h2, L2_psi2, L2_phis2, L2_phib2 = [],[],[],[]

    for t2 in t_steps2:
	    # numerical results
        tt2   = format(t2,'.4f')
        fname2 = tt2+'.npy'
        x_n2, h_n2, psi_n2, phis_n2, phib_n2 = load_binary(data_path2, fname2)

	    # exact standing wave solutions
        h_ex2 = H0 + np.cos(k*x_n2) * ( Aa*np.cos(w*t2) + Bb*np.sin(w*t2) )
        phis_ex2 = np.cos(k*x_n2) * (np.exp(k*H0) + np.exp(-k*H0)) * (-Aa*np.sin(w*t2) + Bb*np.cos(w*t2))/w
        phib_ex2 = np.cos(k*x_n2) * (np.exp(k*0) + np.exp(-k*0)) * (-Aa*np.sin(w*t2) + Bb*np.cos(w*t2))/w

        L2_h2.append(norm(h_ex2-h_n2)) 
        L2_psi2.append(norm(phis_ex2-psi_n2))
        L2_phis2.append(norm(phis_ex2-phis_n2))
        L2_phib2.append(norm(phib_ex2-phib_n2))

fig, (ax1, ax2) = plt.subplots(2,figsize=(9,9),constrained_layout=True)

ax1.set_xlabel('$t$ [s]',fontsize=14)
ax1.set_ylabel('$L^2$ error: $h$',fontsize=14)
ax1.plot(t_steps,L2_h,'ko',label='$h$'+label,markersize=5)
if comparison:
	ax1.plot(t_steps2,L2_h2,'y--',label='$h$'+label2,markersize=3)
ax1.legend(loc='best',fontsize=14)
ax1.grid()

ax2.set_xlabel('$t$ [s]',fontsize=14)
ax2.set_ylabel(r'$L^2$ error: $\phi$',fontsize=14)
ax2.plot(t_steps,L2_psi,'ro',label=r'$\psi$'+label,markersize=5)
ax2.plot(t_steps,L2_phis,'bo',label=r'$\phi|_{z=H_0}$'+label,markersize=3)
ax2.plot(t_steps,L2_phib,'go',label=r'$\phi|_{z=0}$'+label,markersize=5)
if comparison:
    ax2.plot(t_steps2,L2_psi2,'r--',label=r'$\psi$'+label2,markersize=5)
    ax2.plot(t_steps2,L2_phis2,'b-.',label=r'$\phi|_{z=H_0}$'+label2,markersize=3)
    ax2.plot(t_steps2,L2_phib2,'g--',label=r'$\phi|_{z=0}$'+label2,markersize=5)
ax2.legend(loc='best',fontsize=14)
ax2.grid()

if save_figure:
    plt.savefig(save_path,dpi=300)
else:
    plt.show()
