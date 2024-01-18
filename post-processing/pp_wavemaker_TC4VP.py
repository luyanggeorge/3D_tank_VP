import numpy as np
import matplotlib.pyplot as plt
import os.path

save_figure=True

file_exp1='202002/PistonMotion.dat'
file_exp2='202002/PistonVelocity.dat'

data_path='data_VP/MMP/TC4_2D(1ele,phiGLL)_5Jan'
file = os.path.join(data_path,'probes.csv')

label1='Measurements'
label2='Interpolations (MMP)'

if save_figure:
    path=os.path.join(data_path,'figures')
    try: 
        os.makedirs(path, exist_ok = True) 
    except OSError as error:
        print(error)
    figure_name='Wavemaker_MMP.png'
    save_path=os.path.join(path, figure_name)

with open(file_exp1,'r') as fe1:
	t1_exp, R_exp = np.loadtxt(fe1, usecols=(0,1), unpack=True)

with open(file_exp2,'r') as fe2:
	t2_exp, Rt_exp = np.loadtxt(fe2, usecols=(0,1), unpack=True)

with open(file,'r') as fn:
    t_num, R_num, Rt_num = np.loadtxt(fn, usecols=(0,7,8), unpack=True)

plt.figure(num=1, figsize=(12,8),constrained_layout=True)
ax1 = plt.subplot(211)
ax2 = plt.subplot(212)

ax1.set_title('Wavemaker motion', fontsize=20)
ax1.set_xlabel('Time [s]',fontsize=14)
ax1.set_ylabel('R(t) [m]',fontsize=14)
ax1.plot(t1_exp,R_exp,'b-', label=label1)
ax1.plot(t_num,R_num,'r--', label=label2)
ax1.legend(loc='upper right',fontsize=13)
ax1.set_xlim([0,120])
ax1.grid()

ax2.set_title('Wavemaker Velocity', fontsize=20)
ax2.set_xlabel('Time [s]',fontsize=14)
ax2.set_ylabel('dR(t)/dt [m/s]',fontsize=14)
ax2.plot(t2_exp,Rt_exp,'b-',label=label1)
ax2.plot(t_num,Rt_num,'r--',label=label2)
ax2.legend(loc='upper right',fontsize=13)
ax2.set_xlim([0,120])
ax2.grid()

if save_figure:
    plt.savefig(save_path,dpi=300)
else:
    plt.show()