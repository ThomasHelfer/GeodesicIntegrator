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


# ================= Plotting ========================

print("Plotting geodesics ")
N1 = len(data)
cmap = plt.get_cmap('autumn')
colors = [cmap(i) for i in np.linspace(1, 0, N1)]

plt.figure(figsize=(14, 14), dpi=200)
ind = 0
for pos in data:
    plt.plot(pos[1],pos[4],color = colors[ind])
    ind += 1
plt.xlabel(r'$x~[M^{-1}]$')
plt.ylabel(r'$t~[M^{-1}]$')
plt.xlim([0,10])
plt.ylim([0,5])
plt.legend()
plt.grid()
plt.savefig("gedodesic.png",bbox_inches = 'tight')
plt.close()

