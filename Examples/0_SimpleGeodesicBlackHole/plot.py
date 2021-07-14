#########################
# Preamble
########################
import matplotlib
matplotlib.use('agg')

import numpy as np
import pylab as pl
import matplotlib as mpl
import matplotlib.pyplot as plt
import string
import matplotlib.collections as collections
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from scipy.optimize import curve_fit
import os

# ============================
#  Doddy defintions
# ============================

F1 = 30 # Axes font size
F2 = 20 # Legend font size
line = 2.5 # Line width
alp = 1 # alpha
scale= 'linear'
colorshift=3
ExtractionR = 400

# Setting the font structure

rc = mpl.rcParams # Font structure is called rc now
rc['text.usetex'] = True # Tex fonts
rc['font.family'] = 'sans-serif'
rc['font.serif'].insert(0,'cm') # Default font is computer modern for latex
rc['font.size'] = F1
rc['xtick.labelsize'] = 'small'
rc['ytick.labelsize'] = 'small'
rc['legend.fontsize'] = F2

# ==============================
# Read - in 
# ==============================

i = 0 
data = []
filename = "xpos" + str(i).zfill(3)+".csv"
print(filename)
while(os.path.isfile(filename)):
    data.append(np.loadtxt(filename,unpack=True))
    i = i + 1 
    filename = "xpos" + str(i).zfill(3)+".csv"

# ===============================
# Defs 
# ===============================

mass = 1; 
theta = np.linspace(0,2*np.pi,100);
x_horizon = 2*mass*np.sin(theta)
y_horizon = 2*mass*np.cos(theta) 
x_ISCO = 6*mass*np.sin(theta)
y_ISCO = 6*mass*np.cos(theta)
x_Photon = 3*mass*np.sin(theta)
y_Photon = 3*mass*np.cos(theta)

# ================= Plotting ========================

print("Plotting geodesics ")
N1 = len(data)
cmap = plt.get_cmap('autumn')
colors = [cmap(i) for i in np.linspace(1, 0, N1)]

plt.figure(figsize=(14, 14), dpi=200)
plt.plot(x_horizon,y_horizon,label = 'horizon',color = 'black')
plt.plot(x_ISCO,y_ISCO,label = 'ISCO')
plt.plot(x_Photon,y_Photon,label = 'Photon Sphere')
ind = 0 
for pos in data:
    plt.plot(pos[0],pos[1],color = colors[ind])
    ind += 1 
plt.xlabel(r'$x~[M^{-1}]$')
plt.ylabel(r'$y~[M^{-1}]$')
plt.xlim([-30,30])
plt.ylim([-30,30])
plt.legend()
plt.grid()
plt.savefig("gedodesic.png",bbox_inches = 'tight')
plt.close()

refNorm = pos[4][0];

plt.figure(figsize=(14, 14), dpi=200)
for pos in data:
    ref = pos[4][0]
    plt.plot(pos[3],pos[4])#abs(pos[4]-ref)/ref*100)
plt.xlabel(r'$t~[M^{-1}]$')
plt.ylabel(r'$||\partial x/ \partial \tau ||$')
#plt.xlim(-40,max(time_2)/minit[1])
#plt.xlim([-10,10])
#plt.ylim([-1,1])
plt.legend()
plt.grid()
plt.savefig("time_xpos.png",bbox_inches = 'tight')
plt.close()

# Get final angle

plt.figure(figsize=(14, 14), dpi=200)
b_arr = []
for pos in data:
    xdiff = np.diff(pos[0])
    ydiff = np.diff(pos[1])
    theta = -np.arctan2(ydiff[-3],xdiff[-3]) + np.pi
    b = np.abs(pos[1][0])
    b_arr.append(b)
    plt.scatter(b,theta)
b = np.linspace(2.,max(b_arr),100)
plt.plot(b,4./b,color = 'red',label = 'Newtonian')
#plt.ylim([0,np.pi])
plt.xlabel(r'Impact parameter $b$')
plt.ylabel(r'Deflection Angle $\delta\theta$')
plt.savefig('deflection_angle.png')
plt.close()
