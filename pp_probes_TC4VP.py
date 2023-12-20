import numpy as np
import matplotlib.pyplot as plt
import os.path

measurements=[f'202002/WAVE_{i}.dat' for i in range(1,7)]

H0=1
label0='Exprimental data'

save_figure=False

data_path='data_VP/MMP/TC4_2D(one-ele,phiGLL)'
file = os.path.join(data_path,'probes.csv')
label='MMP'


with open(file,'r') as fn:
    tn, p1n, p2n, p3n, p4n, p5n, p6n = np.loadtxt(fn, usecols=(0,1,2,3,4,5,6), unpack=True)

pn_h=[p1n, p2n, p3n, p4n, p5n, p6n]
pn=[np.subtract(pn_h[i],H0) for i in range(len(pn_h))]

te = np.zeros((len(tn),6))
pe = np.zeros((len(tn),6))

for i in range(6):
	with open(measurements[i],'r') as fe:
		te_full, pe_full = np.loadtxt(fe, usecols=(0,1), unpack=True)
	te[:,i]=te_full[:len(tn)]
	pe[:,i]=pe_full[:len(tn)]

if save_figure:
    path=os.path.join(data_path,'figures')
    try: 
        os.makedirs(path, exist_ok = True) 
    except OSError as error:
        print(error)
    figure_name1='Probes_123.png'
    figure_name2='Probes_456.png'
    save_path1=os.path.join(path, figure_name1)
    save_path2=os.path.join(path, figure_name2)


plt.figure(num=1, figsize=(13,9),constrained_layout=True)
ax1 = plt.subplot(311)
ax2 = plt.subplot(312)
ax3 = plt.subplot(313)

ax1.set_title('Probe 1: x=10m', fontsize=20)
ax1.set_xlabel('Time [s]',fontsize=14)
ax1.set_ylabel('$\eta$ [m]',fontsize=14)
ax1.plot(te[:,0],pe[:,0],'r-', label=label0)
ax1.plot(tn,pn[0],'b--', label=label)
#ax1.set_xlim([0,120])
ax1.legend(loc='upper left',fontsize=13)
ax1.grid()

ax2.set_title('Probe 2: x=20m', fontsize=20)
ax2.set_xlabel('Time [s]',fontsize=14)
ax2.set_ylabel('$\eta$ [m]',fontsize=14)
ax2.plot(te[:,1],pe[:,1],'r-',label=label0)
ax2.plot(tn,pn[1],'b--',label=label)
#ax2.set_xlim([0,120])
ax2.legend(loc='upper left',fontsize=13)
ax2.grid()

ax3.set_title('Probe 3: x=40m', fontsize=20)
ax3.set_xlabel('Time [s]',fontsize=14)
ax3.set_ylabel('$\eta$ [m]',fontsize=14)
ax3.plot(te[:,2],pe[:,2],'r-',label=label0)
ax3.plot(tn,pn[2],'b--',label=label)
#ax3.set_xlim([0,120])
ax3.legend(loc='upper left',fontsize=13)
ax3.grid()

if save_figure:
    plt.savefig(save_path1,dpi=300)
else:
    plt.show()

plt.figure(num=2, figsize=(13,9),constrained_layout=True)
ax1 = plt.subplot(311)
ax2 = plt.subplot(312)
ax3 = plt.subplot(313)

ax1.set_title('Probe 4: x=49.5m', fontsize=20)
ax1.set_xlabel('Time [s]',fontsize=14)
ax1.set_ylabel('$\eta$ [m]',fontsize=14)
ax1.plot(te[:,3],pe[:,3],'r-', label=label0)
ax1.plot(tn,pn[3],'b--', label=label)
#ax1.set_xlim([0,120])
ax1.legend(loc='upper left',fontsize=13)
ax1.grid()

ax2.set_title('Probe 5: x=50m', fontsize=20)
ax2.set_xlabel('Time [s]',fontsize=14)
ax2.set_ylabel('$\eta$ [m]',fontsize=14)
ax2.plot(te[:,4],pe[:,4],'r-',label=label0)
ax2.plot(tn,pn[4],'b--',label=label)
#ax2.set_xlim([0,120])
ax2.legend(loc='upper left',fontsize=13)
ax2.grid()

ax3.set_title('Probe 6: x=54m', fontsize=20)
ax3.set_xlabel('Time [s]',fontsize=14)
ax3.set_ylabel('$\eta$ [m]',fontsize=14)
ax3.plot(te[:,5],pe[:,5],'r-',label=label0)
ax3.plot(tn,pn[5],'b--',label=label)
#ax3.set_xlim([0,120])
ax3.legend(loc='upper left',fontsize=13)
ax3.grid()

if save_figure:
    plt.savefig(save_path2,dpi=300)
else:
    plt.show()