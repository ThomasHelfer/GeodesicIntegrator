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

# Set a color stucture

# ===============================
# Defs 
# ===============================


pos = np.loadtxt("xpos.csv",unpack=True)

print(pos)
mass = 1; 
theta = np.linspace(0,2*np.pi,100);
x_horizon = 2*mass*np.sin(theta)
y_horizon = 2*mass*np.cos(theta) 
x_ISCO = 6*mass*np.sin(theta)
y_ISCO = 6*mass*np.cos(theta)

plt.figure(figsize=(14, 14), dpi=200)
plt.plot(pos[0],pos[1],label = "Geodesic" )
plt.plot(x_horizon,y_horizon,label = 'horizon',color = 'black')
plt.plot(x_ISCO,y_ISCO,label = 'ISCO')
plt.xlabel(r'$x~[M^{-1}]$')
plt.ylabel(r'$y~[M^{-1}]$')
#plt.xlim(-40,max(time_2)/minit[1])
#plt.xlim([-10,10])
#plt.ylim([-10,10])
plt.legend()
plt.grid()
plt.savefig("gedodesic.png",bbox_inches = 'tight')
plt.close()


plt.figure(figsize=(14, 14), dpi=200)
plt.plot(pos[3],pos[1],label = "x position " )
#plt.xlabel(r'$x/m_{init}$')
#plt.ylabel(r'$r\psi_4 m_{init}$')
#plt.xlim(-40,max(time_2)/minit[1])
#plt.xlim([-10,10])
#plt.ylim([-10,10])
plt.legend()
plt.grid()
plt.savefig("time_xpos.png",bbox_inches = 'tight')
plt.close()
