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


fig = plt.figure(figsize=(1.5*6.3,1.5*2.7))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

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

dat = np.genfromtxt("data/qtrans/k4_3000ens.txt")
q, corner, middle = np.hsplit(dat, 3)

ax1.plot(q, corner, "o", color=martaGold, ms=5.0, alpha=0.9)
ax1.plot(q, middle, "o", color=martaPurple, ms=5.0, alpha=0.9)
ax2.plot(q, corner-middle, "o", ms=5.0, color=martaRed)
ax2.axhline(0.0, color="k")
plt.show() 
