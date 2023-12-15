# Post-processing code for TC1 using the VP approach
# Validate the numerical results against the standing wave solutions of linear potential-flow eqs
# Making an animation flame by flame

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
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

slice_y = 0

data_path ='data_VP/MMP/TC1_test4_2D(one-ele,phiGLL)'
gif_name='TC1_test4(one_ele,hatphi_GLL).gif'

interval=1
#------------------------#

# parameters
k = 2*np.pi*m/Lx
g  = 9.81
kH0 = np.arccosh(0.5*g)
H0 = kH0/k
w = np.sqrt(2*k*np.sinh(k*H0))
amp_eta=np.sqrt(Aa*Aa + Bb*Bb)
amp_phi=amp_eta*(np.exp(k*H0)+np.exp(-k*H0))/w

path=os.path.join(data_path,'figures')
try: 
    os.makedirs(path, exist_ok = True)
except OSError as error: 
    print(error) 

# get t from the energy file
energy_file = os.path.join(data_path,'energy.csv')
with open(energy_file,'r') as e_f:
    t_steps_full = np.loadtxt(e_f, usecols=0)
t_steps = t_steps_full[::interval]

fig, (ax1, ax2) = plt.subplots(2,figsize=(7,7),constrained_layout=True)

line_he, = ax1.plot([], [], 'k-',label='Exact')
line_hn, = ax1.plot([], [], 'r--',label='Numerical')
line_phise, = ax2.plot([], [], 'k-',label='$z=H_0$, exact')
line_phisn1,= ax2.plot([], [], 'r--',label='$z=H_0$, $\psi$')
line_phisn2,= ax2.plot([], [], 'b-.',label='$z=H_0$, $\phi$')
line_phibe, = ax2.plot([], [], 'g-',label='$z=0$, exact')
line_phibn, = ax2.plot([], [], 'c-.',label='$z=0$, $\phi$')
time_text=ax1.text( 0.8, 0.9, '', transform=ax1.transAxes, fontsize='large',
                   bbox=dict(boxstyle="round", ec=(1., 0.5, 0.5), fc=(1., 0.8, 0.8)) )
ax1.grid()
ax2.grid()

def init():
    ax1.set_xlim(0, Lx)
    ax1.set_ylim(-1.1*amp_eta+H0, 1.1*amp_eta+H0)
    ax1.set_title(f'$h(x,y={slice_y},t)$',fontsize='x-large')
    ax1.set_xlabel('$x$',fontsize='large')
    ax1.set_ylabel('$h$',fontsize='large')
    ax1.legend(loc='upper left',fontsize='large')
    
    ax2.set_xlim(0, Lx)
    ax2.set_ylim(-1.1*amp_phi, 1.1*amp_phi)
    ax2.set_title(r'$\phi(x,y=0,z,t)$',fontsize='x-large')
    ax2.set_xlabel('$x$',fontsize='large')
    ax2.set_ylabel(r'$\phi$',fontsize='large')
    ax2.legend(loc='upper left',fontsize='large')
    return line_he,line_hn,line_phise,line_phisn1,line_phisn2,line_phibe,line_phibn,

def animate(t):
    # numerical results
    tt   = format(t,'.4f')
    fname = tt+'.npy'
    x_n, h_n, psi_n, phis_n, phib_n = load_binary(data_path, fname)

    # exact standing wave solutions
    h_ex = H0 + np.cos(k*x_n) * ( Aa*np.cos(w*t) + Bb*np.sin(w*t) )
    phis_ex = np.cos(k*x_n) * (np.exp(k*H0) + np.exp(-k*H0)) * (-Aa*np.sin(w*t) + Bb*np.cos(w*t))/w
    phib_ex = np.cos(k*x_n) * (np.exp(k*0) + np.exp(-k*0)) * (-Aa*np.sin(w*t) + Bb*np.cos(w*t))/w

    # plot
    line_he.set_data(x_n, h_ex)
    line_hn.set_data(x_n, h_n)
    line_phise.set_data(x_n, phis_ex)
    line_phisn1.set_data(x_n, psi_n)
    line_phisn2.set_data(x_n, phis_n)
    line_phibe.set_data(x_n, phib_ex)
    line_phibn.set_data(x_n, phib_n)

    time_text.set_text('$t = %s$' %tt)

    return line_he,line_hn,line_phise,line_phisn1,line_phisn2,line_phibe,line_phibn,

anim = FuncAnimation(fig, animate, init_func=init, frames=t_steps, interval=200, blit=True)

anim.save(os.path.join(path, gif_name), writer='pillow',dpi=200)

plt.close()