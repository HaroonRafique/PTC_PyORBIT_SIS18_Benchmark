#-----------------------------------------------------------------------
# This plotting script is used to calculate the single particle tune
# using the particle co-ordinates which are recorded once per turn
# Note that this requires access to the NAFF libraries, kindly made 
# available by Nikolaos Karastathis (CERN BE-ABP-HSI)
# The data are plotted crudely on top of that provided by Giuliano 
# Franchetti here: 
# https://web-docs.gsi.de/~giuliano/research_activity/trapping_benchmarking/main.html
#-----------------------------------------------------------------------

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

# PARAMETERS
horizontal = 0    
n_part = 100
n_turns = 1024
min_particle = 1

intensity = 'I=2.95e10'
step = '4'
plane =''

if horizontal:
    plane = 'x'
else:
    plane = 'y'
    
extra = ''
case = 'Step' + step + '_' + plane + '_' + intensity

if len(extra) > 0:
    case = case + '_' + extra

Qx = 4.338
Qy = 3.2

if '4' in step:
    Qx = 4.3504
else:
    Qx = 4.338

# NORMALISATION FOR SIS18
beta = 0.15448
gamma = 1.012149995

bety0=13.47433765
alfy0=0.4264503497
gamy0 = (1+alfy0*alfy0)/bety0
eg_y = (9.3e-6)/4
sig_y = np.sqrt(eg_y * bety0)

betx0=12.79426135
alfx0=1.283306757 
gamx0 = (1+alfx0*alfx0)/betx0
eg_x =  (12.57e-6)/4
sig_x = np.sqrt(eg_x * betx0)

const_x = np.sqrt(betx0 * gamx0) 
const_y = np.sqrt(bety0 * gamy0) 

const =1.
if horizontal:
    const=const_x/sig_x
else:
    const=const_y/sig_y

# IMPORT DATA
filename = 'mainbunch'

files = glob.glob(filename + '*.mat')
files.sort()
df_data=pd.DataFrame({})

for i, file in enumerate(files):
    particles = sio.loadmat(file)
    data=(zip(particles['particles'][0][0]['x'].flatten().tolist(),particles['particles'][0][0]['xp'].flatten().tolist(),particles['particles'][0][0]['y'].flatten().tolist(),particles['particles'][0][0]['yp'].flatten().tolist(),particles['particles'][0][0]['z'].flatten().tolist(),particles['particles'][0][0]['dE'].flatten().tolist(),particles['particles'][0][0]['ParticleIdNumber'].flatten().tolist(),particles['particles'][0][0]['macrosize'].flatten().tolist()))
    data=np.array(sorted(data, key=lambda x:x[-1]))
    df=pd.DataFrame(np.array(data)[:,:7],index=np.array(data).astype('int64')[:,6], columns=['x','xp','y','yp','z','dE','macrosize'])
    df_data=df_data.append(df)

if min_particle%2:
    turns_half=((n_turns/2)-(min_particle+1))
else:
    turns_half=(n_turns/2-(min_particle))

# CALCULATE TUNES USING NAFF
tunes=[]
for i in xrange(min_particle,n_part):
    if horizontal:
        tunes.append(naff(data=np.array(df_data['x'].iloc[i::n_part]), turns=turns_half, nterms=1)[0][1])
    else:
        tunes.append(naff(data=np.array(df_data['y'].iloc[i::n_part]), turns=turns_half, nterms=1)[0][1])

max_turn=1023
z=[]

# CALCULATE X-AXIS DATA
for i in xrange(min_particle,n_part):
    if horizontal:
        z.append( const * np.array(df_data['x'].iloc[i:(i+1):n_part])[0]) 
    else:
        z.append( const * np.array(df_data['y'].iloc[i:(i+1):n_part])[0]) 

print 'length of x-axis data = ',len(z)
print 'length of tunes = ',len(tunes)

fig, ax = plt.subplots()

plot_title=''
if horizontal:
	plot_title = 'SIS18_'+ plane +'_Distn Qx = ' + str(Qx) + '_' + case
	ax.set(xlabel='x [sigma_x]', ylabel='Qx', title=plot_title)
	plot_name = 'SIS18_' + case + '_Qx.png'
	img = plt.imread("step4_x.png")
	# ~ ax.imshow(img, zorder=0, origin='upper', aspect='auto', extent=[0.0, 5.0, 0.2325, 0.3375]) #Step 2
	# ~ ax.imshow(img, zorder=0, origin='upper', aspect='auto', extent=[0.0, 7.0, 0.2325, 0.345]) #Step 3
	ax.imshow(img, zorder=0, origin='upper', aspect='auto', extent=[0.0, 6.0, 0.2325, 0.3505]) #Step 4

else:
	plot_title = 'SIS18_'+ plane +'_Distn Qy = ' + str(Qy) + '_' + case
	ax.set(xlabel='y [sigma_y]', ylabel='Qy', title=plot_title)
	plot_name = 'SIS18_' + case + '_Qy.png'
	ax.set_xlim(0,5)
	ax.set_ylim(0.05, 0.20625)		
	img = plt.imread("step4_y.png")
	# ~ ax.imshow(img, zorder=0, origin='upper', aspect='auto', extent=[0.0, 5.0, 0.05, 0.205]) #Step 2
	# ~ ax.imshow(img, zorder=0, origin='upper', aspect='auto', extent=[0.0, 5.0, 0.05, 0.20625]) #Step 3
	ax.imshow(img, zorder=0, origin='upper', aspect='auto', extent=[0.0, 4.62, 0.05, 0.205]) #Step 4

ax.grid()
ax.scatter(z, tunes, marker='.', zorder=1, color='m', label='PTC PyORBIT')
ax.legend()

fig.savefig(plot_name, dpi=600)

# Save file
outfile_name = 'SIS18_' + case + '_Tunes_OnData.txt'
np.savetxt(outfile_name, zip(z,tunes), fmt="%3.6e %3.6e")
print 'Tune plot completed'
