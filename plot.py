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

plt.figure(figsize=(14, 14), dpi=200)
plt.plot(pos[0],pos[1],label = "x position " )
#plt.xlabel(r'$x/m_{init}$')
#plt.ylabel(r'$r\psi_4 m_{init}$')
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
