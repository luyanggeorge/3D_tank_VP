# Post-processing code for TC1 using the VP approach or the old approach (flag: auto_WF)
import numpy as np
import matplotlib.pyplot as plt
import os.path

#------ User Input ------#
save_figure=True
figure_name='TC1_energy_comparison_SE(dt-dt2).png'

rel_err=True
auto_WF=False

data_path ='data_Re/SE/TC1_2D_test_dt'
label='$\Delta t$'

comparison=True
data_path2='data_Re/SE/TC1_2D_test_dt2'
label2='$\Delta t/2$'
#------------------------#

file = os.path.join(data_path,'energy.csv')
with open(file,'r') as f:
    if auto_WF:
        time, Ek, Ep, E_tot = np.loadtxt(f, usecols=(0,1,2,3), unpack=True)
    else:
        time, E_tot = np.loadtxt(f, usecols=(0,1), unpack=True)

if rel_err:
    E_err = np.array([(E_tot[i]-E_tot[0])/E_tot[0] for i in range(len(E_tot))])

if save_figure:
    path=os.path.join(data_path,'figures')
    try: 
        os.makedirs(path, exist_ok = True) 
    except OSError as error:
        print(error)
    save_path=os.path.join(path, figure_name)

if comparison:
    file2 = os.path.join(data_path2,'energy.csv')
    with open(file2,'r') as f2:
        if auto_WF:
            time2, Ek2, Ep2, E_tot2 = np.loadtxt(f2, usecols=(0,1,2,3), unpack=True)
        else:
            time2, E_tot2 = np.loadtxt(f2, usecols=(0,1), unpack=True)

    if rel_err:
        E_err2 = np.array([(E_tot2[i]-E_tot2[0])/E_tot2[0] for i in range(len(E_tot2))])

fig, ax1=plt.subplots()
fig.set_size_inches(8,4)
fig.set_tight_layout(True)

ax1.set_title('Energy variations', fontsize=18)
ax1.set_xlabel('$t$ [s]',fontsize=14)

#ax1.ticklabel_format(axis='y',scilimits=(0,0))
#ax1.plot(time,Ek,'bo',label='Kinetic energy')
#ax1.plot(time,Ep,'bo',label='Potential energy')
if rel_err:
    ax1.set_ylabel('$[E(t)-E(t_0)]/E(t_0)$',fontsize=14)
    ax1.plot(time,E_err,'ko',label=label,markersize=5)
else:
    ax1.set_ylabel('$E(t)$',fontsize=14)
    ax1.plot(time,E_tot,'ko',label=label,markersize=5)


if comparison:
    if rel_err:
        ax1.plot(time2,E_err2,'ro',label=label2,markersize=5)
        ax1.plot(time2,E_err2*2,'yo',label=r'$\Delta t/2 \times 2$',markersize=3)
    else:
        ax1.plot(time2,E_tot2,'ro',label=label2,markersize=5)
    ax1.legend(loc='upper left',fontsize=14)
    
ax1.grid()

if save_figure:
    plt.savefig(save_path,dpi=300)
else:
    plt.show()