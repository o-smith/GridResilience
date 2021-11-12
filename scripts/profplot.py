import glob, os
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt 
from matplotlib import colors
# from ternary.helpers import simplex_iterator
# import ternary 
from scipy.optimize import curve_fit
from scipy.ndimage.filters import gaussian_filter

#Set up plot style
font = {'size'   : 20}
matplotlib.rc('font', **font)
matplotlib.rc('font', serif='Computer Modern Roman')


#Define colours 
martaRed = "#c24c51" 
martaGreen = "#54a666"
martaBlue = "#4c70b0"
martaPurple = "#7f70b0"
martaGold = "#ccb873"
Icolour = '#DB2420' #Central line red
Ocolour = '#00A0E2' #Victoria line blue
O2colour = '#868F98' #Jubilee line grey
D2colour = '#F386A0' #Hammersmith line pink 
D4colour = '#97015E' #Metropolitan line magenta
D6colour = '#B05F0F' #Bakerloo line brown
D5colour = '#00843D' #District line green

dat = np.genfromtxt("data/new/q1_15_45_0_time.txt") 
x, y, z, z2 = np.hsplit(dat, 4) 
x = np.array(x)[:,0] 
z = np.array(z)[:,0]
print(np.shape(x))

fig = plt.figure()
ax = fig.add_subplot(111)


ax.plot(x,z, color=martaRed, lw=3.0) 
ax.axvline(0.94+0.1, color="k", lw=3.0, alpha=0.7)
ax.fill_between(x, z, z*0.0, color=martaRed, alpha=0.4)
ax.set_ylim([0,30])
ax.set_xlim([0.1,2.5])

plt.tight_layout()
plt.show() 
















