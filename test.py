import functions as fn
import modules as m
import numpy as np
from matplotlib import pyplot as plt
import sys

#Define some a priori known parameters (in practice, these variables will change according to the pair)
rows = 6752
cols =  1700
master = 'test_data/SAO1A_20190820_HH'
slave = 'test_data/SAO1A_20191124_HH'

#Read the imput sample data
x = np.fromfile('test_data/X_Regr.dat', dtype= np.float64).reshape(rows, cols)
y = np.fromfile('test_data/Y_Regr.dat', dtype= np.float64).reshape(rows, cols)
z = np.fromfile('test_data/Z_Regr.dat', dtype= np.float64).reshape(rows, cols)

coh = np.fromfile('test_data/coh.dat', dtype= np.float32)
coh = coh.reshape(rows, cols)

fas = np.fromfile('test_data/fr.dat',  dtype= np.float32)
fas = fas.reshape(rows, cols)

m1, parm1 = fn.read_params_gmt(master)
s1, pars1 = fn.read_params_gmt(slave)

#Run the Orbital Corrector
correct = m.mainCorrector(x,y,z,coh,fas,m1,s1,4,4, ramp = True).runCorrector()

#Plot the results

ifg = '20190820_20191124'

fig, axs = plt.subplots(nrows=1, ncols=2, sharex=True)
fig.suptitle(f'Results for pair {ifg}')
a = axs[0].imshow(fas, cmap = 'hsv_r')
axs[0].set_title("Original")
clb = fig.colorbar(a,ticks=[-3,  3], ax = axs[0])
clb.ax.set_yticklabels([r'$-\pi$',r'$\pi$'])
clb.ax.tick_params(labelsize=8,size=0)
b = axs[1].imshow(correct[0], cmap = 'hsv_r')
axs[1].set_title("Corrected")
clb = fig.colorbar(b,ticks=[-3,  3], ax = axs[1])
clb.ax.set_yticklabels([r'$-\pi$',r'$\pi$'])
clb.ax.tick_params(labelsize=8,size=0)
plt.show()

