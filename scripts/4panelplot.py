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


fig = plt.figure(figsize=(8,6.5)) 
ax1 = fig.add_subplot(222)
ax2 = fig.add_subplot(224)

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

dat = np.genfromtxt("data/kcont/n64_q0.txt")
k, rho = np.hsplit(dat, 2) 
ax1.plot(k, rho, "-", lw=4.0, color=martaGold, label="$n=64$")
ks = k

dat = np.genfromtxt("data/kcont/n32_q0.txt")
k, rho = np.hsplit(dat, 2) 
ax1.plot(k, rho, "-", lw=4.0, color=martaRed, label="$n=32$")

# dat = np.genfromtxt("data/kcont/n24_q0.txt")
# k, rho = np.hsplit(dat, 2) 
# ax1.plot(k, rho, "-", lw=3.0, color=martaBlue)

dat = np.genfromtxt("data/kcont/n16_q0.txt")
k, rho = np.hsplit(dat, 2) 
ax1.plot(k, rho, "-", lw=4.0, color=martaGreen, label="$n=16$")

ax1.plot(ks, ks/(ks-1), "k--", lw=3, alpha=0.6)
ax1.set_yticks([1,2])
ax1.set_xticks([0,60])
ax1.set_ylabel("$\\rho$", rotation=0, fontsize=22)
ax1.set_xlabel("$K$", fontsize=20, labelpad=-15)
ax1.legend(prop={'size': 12})

# dat = np.genfromtxt("data/qtrans/k4_2500ens.txt")
# q, acc, acm = np.hsplit(dat, 3)
# plt.plot(q, acc-acm, ".", color=martaRed, ms=9.0)
# plt.xlim([0,1])
# plt.xticks([0,1])
# plt.ylim([-0.4,0.3])
# plt.yticks([-0.4, 0.0, 0.3])
# plt.axhline(color="k")
# plt.xlabel("$q$", fontsize=24)
# plt.ylabel("$\\phi$", rotation=0, labelpad=20)
# plt.imshow(dat)
# plt.imshow(dat, cmap="Blues", aspect="auto", extent=[0.01,3,1,49], vmin=0, vmax=1)
# ax2.set_yticks([1,49])
# ax2.set_xticks([0,3])
# ax2.set_ylabel("$d$", rotation=0, fontsize=22, labelpad=-8)
# ax2.set_xlabel("$\\alpha/\\alpha_{\\ast}$", fontsize=22, labelpad=-15)
# cax = fig.add_axes([0.27, 0.8, 0.5, 0.05])
# plt.colorbar(shrink=1)
plt.tight_layout() 
plt.show() 




# plt.subplots_adjust(left=0.2, wspace=0.4, hspace=0.6)
# # plt.tight_layout() 
# plt.gcf().subplots_adjust(bottom=0.15)
