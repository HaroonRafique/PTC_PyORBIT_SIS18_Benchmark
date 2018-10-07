import os
import sys
import numpy as np
import scipy.io as sio
import glob
import pandas as pd
from matplotlib import cm
from scipy.interpolate import griddata
import numpy.ma as ma
import pylab as plt
import sys
sys.path.append("/eos/user/n/nkarast/SWAN/PyNAFF/")
from PyNAFF import naff
import matplotlib.pyplot as plt

file_in = 'output.mat'
main_label = 'Long Distn Test: Voltage'

particles = dict()
sio.loadmat(file_in, mdict=particles)

fig1, ax1 = plt.subplots()

img = plt.imread("step9.png")
# Synch. Osc.
ax1.imshow(img, zorder=0, origin='upper', aspect='auto', extent=[0.0, 100., 1.0, 2.2]) #Step 9
ax1.plot( ((particles['turn'][0])/(1000.)), (particles['epsn_x'][0]/particles['epsn_x'][0][0]), label='PTC PyORBIT', linewidth = 2);
ax1.set_xlabel('Synch. Osc. [-]');
ax1.set_xlim(0, 100)

ax1.set_ylabel('eps_x / eps_x_0 [-]')
ax1.set_ylim(1, 2.2)		

ax1.set_title('SIS18 Step 9');
ax1.grid(True);
ax1.legend()

figname = 'SIS18_Step9.png'	
fig1.savefig(figname);	
