# Compute the rate of convergence using Atiken Extrapolation
# Using 5x5 points

import numpy as np
from numpy.linalg import norm
from matplotlib import pyplot as plt
import os.path

# -----------------------
g = 9.81
H0 = 1.0
lamb = 2.0                      # Wavelength
k = 2*np.pi/lamb                # Wave number
w = np.sqrt(g*k*np.tanh(k*H0))  # Wave frequency
Tw = 2*np.pi/w

n_z  = '8'
nCG  = 'CG1'
ave_t0   = 0.002 # the start time for averaging
ave_t1   = 7

# 1->5 from coarse to fine meshes
# three groups: A{0.2,0.1,0.05}, B{0.1,0.05,0.025}, C{0.05,0.025,0.0125}

if n_z=='4':
	data_path1 = 'data_VP/MMP/TC2_nz4_CG1/h02_dt2'
	data_path2 = 'data_VP/MMP/TC2_nz4_CG1/h01_dt2'
	data_path3 = 'data_VP/MMP/TC2_nz4_CG1/h005_dt2'
	data_path4 = 'data_VP/MMP/TC2_nz4_CG1/h0025_dt2'
	data_path5 = 'data_VP/MMP/TC2_nz4_CG1/h00125_dt2'
elif n_z=='8':
	data_path1 = 'data_VP/MMP/TC2_nz8_CG1/h02_dt2'
	data_path2 = 'data_VP/MMP/TC2_nz8_CG1/h01_dt2'
	data_path3 = 'data_VP/MMP/TC2_nz8_CG1/h005_dt2'
	data_path4 = 'data_VP/MMP/TC2_nz8_CG1/h0025_dt2'
	data_path5 = 'data_VP/MMP/TC2_nz8_CG1/h00125_dt2'

save_figure = True

figure_name = 'Atiken_MMP(vertices,'+nCG+f',nz{n_z}).png'
# -----------------------

if save_figure:
    path=os.path.join(data_path1,'figures')
    try: 
        os.makedirs(path, exist_ok = True) 
    except OSError as error: 
        print(error)
    save_path = os.path.join(path, figure_name)

def load_data(file_path, file_name):
	"""load data from a binary file"""
	binary_file = os.path.join(file_path, file_name)
	with open(binary_file,'rb') as bf:
		d = np.load(bf)
	#return d # full set
	#return d[[0,8,-9,-1]] # four vertices 2*2
	#return d[[0,4,8,180,184,188,-9,-5,-1]] # 3*3
	#return np.concatenate([ d[0:9:2], d[90:99:2], d[180:189:2], d[270:279:2], d[360:369:2] ]) # 5*5
	return np.concatenate([d[i*9:(i*9+9):2] for i in range(0,41,2)]) # 5*21, vertices

def L2norm(array1, array2, weight):
	"""compute the approximated L2-norm of two group of results"""
	diff = (array1-array2)**2
	appro = diff*weight
	integral = np.sum(appro)
	return np.sqrt(integral)

left_right = np.array([0.25, 0.5, 0.5, 0.5, 0.25])
middle = np.array([0.5, 1, 1, 1, 0.5])
weight = np.concatenate([left_right, np.tile(middle,19), left_right])
#weight = np.tile([1],4*20)
print(weight)

t_steps = np.arange(0.002,7.002,0.002)

# rate of convergence 's'
L2_12 = []
L2_23 = []
L2_34 = []
L2_45 = []

Li_12 = []
Li_23 = []
Li_34 = []
Li_45 = []

for t in t_steps:
	# numerical results
	tt   = format(t,'.4f')
	fname = tt+'.npy'

	h1 = load_data(data_path1,fname)
	h2 = load_data(data_path2,fname)
	h3 = load_data(data_path3,fname)
	h4 = load_data(data_path4,fname)
	h5 = load_data(data_path5,fname)
	
	L2_12.append(L2norm(h1,h2,weight))
	L2_23.append(L2norm(h2,h3,weight))
	L2_34.append(L2norm(h3,h4,weight))
	L2_45.append(L2norm(h4,h5,weight))

	Li_12.append(np.max(np.abs(h1-h2)))
	Li_23.append(np.max(np.abs(h2-h3)))
	Li_34.append(np.max(np.abs(h3-h4)))
	Li_45.append(np.max(np.abs(h4-h5)))

# Aitken Extrapolation
s_L2_A = np.log2(np.array(L2_12)/np.array(L2_23))
s_L2_B = np.log2(np.array(L2_23)/np.array(L2_34))
s_L2_C = np.log2(np.array(L2_34)/np.array(L2_45))

s_Li_A = np.log2(np.array(Li_12)/np.array(Li_23))
s_Li_B = np.log2(np.array(Li_23)/np.array(Li_34))
s_Li_C = np.log2(np.array(Li_34)/np.array(Li_45))

# average values:
N_t0 = round(ave_t0/0.002)-1 
N_t1 = round(ave_t1/0.002)-1

ave_s2_A = np.average(s_L2_A[N_t0:N_t1])
ave_s2_B = np.average(s_L2_B[N_t0:N_t1])
ave_s2_C = np.average(s_L2_C[N_t0:N_t1])

ave_si_A = np.average(s_Li_A[N_t0:N_t1])
ave_si_B = np.average(s_Li_B[N_t0:N_t1])
ave_si_C = np.average(s_Li_C[N_t0:N_t1])

print(f'Average values (between t = {t_steps[N_t0]} and {t_steps[N_t1]} [s]):')
print(f's_L2_123 = {ave_s2_A:.4f}; s_Li_123 = {ave_si_A:.4f}')
print(f's_L2_234 = {ave_s2_B:.4f}; s_Li_234 = {ave_si_B:.4f}')
print(f's_L2_345 = {ave_s2_C:.4f}; s_Li_345 = {ave_si_C:.4f}')

# plot
fig, ax = plt.subplots(figsize=(10, 5),constrained_layout=True)
ax.set_title(f'Rate of Convergence ($n_z={n_z}$)', fontsize=18)
ax.set_xlabel('$t$ [s]',fontsize=15)
ax.set_ylabel('$s(t)$',fontsize=15)
ax.set_xlim(0,t_steps[-1])
ax.plot(t_steps, s_L2_A, 'r-', label=r'$L^2, \{h_{0.2}, h_{0.1}, h_{0.05}\}$')
ax.plot(t_steps, s_Li_A, 'r--',label=r'$L^{\infty}, \{h_{0.2}, h_{0.1}, h_{0.05}\}$')
ax.plot(t_steps, s_L2_B, 'g-', label=r'$L^2, \{h_{0.1}, h_{0.05}, h_{0.025}\}$')
ax.plot(t_steps, s_Li_B, 'g--',label=r'$L^{\infty}, \{h_{0.1}, h_{0.05}, h_{0.025}\}$')
ax.plot(t_steps, s_L2_C, 'b-', label=r'$L^2, \{h_{0.05}, h_{0.025}, h_{0.0125}\}$')
ax.plot(t_steps, s_Li_C, 'b--',label=r'$L^{\infty}, \{h_{0.05}, h_{0.025}, h_{0.0125}\}$')
ax.legend(loc='best',fontsize=13, ncol=3)
ax.grid()

if save_figure:
    plt.savefig(save_path, dpi=300)
else:
    plt.show()