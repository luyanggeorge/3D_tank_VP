import numpy as np
import matplotlib.pyplot as plt
import os.path

comparison=True
save_figure=True
rel_err=True

data_path='data_VP/MMP/TC1_test4_2D(one-ele,phiGLL)'
file = os.path.join(data_path,'energy.csv')
label=r'2D(one-ele, $\hat{\phi}=$GLL)'

with open(file,'r') as f:
    time, Ek, Ep, E_tot = np.loadtxt(f, usecols=(0,1,2,3), unpack=True)

if rel_err:
    E_err = np.array([(E_tot[i]-E_tot[0])/E_tot[0] for i in range(len(E_tot))])

if save_figure:
    path=os.path.join(data_path,'figures')
    try: 
        os.makedirs(path, exist_ok = True) 
    except OSError as error:
        print(error)
    figure_name='TC1_energy_comparison(GLL-1).png'
    save_path=os.path.join(path, figure_name)

if comparison:
    data_path2='data_VP/MMP/TC1_test1_dt'
    file2 = os.path.join(data_path2,'energy.csv')
    label2=r'3D(mul-ele, $\hat{\phi}=1$)'

    with open(file2,'r') as f2:
        time2, Ek2, Ep2, E_tot2 = np.loadtxt(f2, usecols=(0,1,2,3), unpack=True)

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
        ax1.plot(time2,E_err2,'ro',label=label2,markersize=4)
        #ax1.plot(time2,E_err2*4,'go',label=r'$Error(\Delta t/2) \times 4$',markersize=4)
    else:
        ax1.plot(time2,E_tot2,'ro',label=label2,markersize=5)
    ax1.legend(loc='upper left',fontsize=14)
    
ax1.grid()

if save_figure:
    plt.savefig(save_path,dpi=300)
else:
    plt.show()