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
font = {'size'   : 18}
matplotlib.rc('font', **font)
matplotlib.rc('font', serif='Computer Modern Roman')


fig = plt.figure(figsize=(6.3,2.7)) 
ax1 = fig.add_subplot(122)


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

dat = np.genfromtxt("data/kcont/n64_q1.txt")
k, rho, _ = np.hsplit(dat, 3) 
ax1.plot(k, rho, "o", color=martaGold, ms=4, label="$n=64$, $q=1$")
ks = k

dat = np.genfromtxt("data/kcont/n32_q1.txt")
k, rho, _ = np.hsplit(dat, 3) 
ax1.plot(k, rho, "o", color=martaRed, ms=4, label="$n=32$, $q=1$")

# dat = np.genfromtxt("data/kcont/n24_q0.txt")
# k, rho = np.hsplit(dat, 2) 
# ax1.plot(k, rho, "-", lw=3.0, color=martaBlue)

# dat = np.genfromtxt("data/kcont/n60_q01.txt")
# k, rho, _ = np.hsplit(dat, 3) 
# ax1.plot(k, rho, "*", color=martaGreen, label="$n=16$")

dat = np.genfromtxt("data/kcont/n70_q1.txt")
k, rho, _ = np.hsplit(dat, 3) 
ax1.plot(k, rho, "o", color=martaGreen, ms=4, label="$n=70$, $q=1$")

dat = np.genfromtxt("data/kcont/n70_q5.txt")
k, rho, _ = np.hsplit(dat, 3) 
ax1.plot(k, rho, "o", color=martaBlue, ms=4, label="$n=16$, $q=0.5$")

dat = np.genfromtxt("data/kcont/n70_q06.txt")
k, rho, _ = np.hsplit(dat, 3) 
ax1.plot(k, rho, "o", color=martaPurple, ms=4, label="$n=16$, $q=0.06$")

ax1.plot(ks, ks/(ks-1), "k--", lw=4.5, alpha=0.8)
ax1.set_yticks([1,1.4])
ax1.set_xticks([0,70])
ax1.set_ylabel("$\\rho$", rotation=0, fontsize=22)
ax1.set_xlabel("$K$", fontsize=20, labelpad=-15)
ax1.legend(prop={'size': 8})


plt.subplots_adjust(left=0.2, wspace=0.7, hspace=0.6)
plt.tight_layout()
plt.show() 