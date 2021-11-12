import numpy as np 
import matplotlib 
from random import random
from numpy.random import rand 
import matplotlib.pyplot as plt 
import matplotlib.colors as mcolors
from scipy.stats import gaussian_kde

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

p1p = np.genfromtxt("data/p1p.txt")
p1g = np.genfromtxt("data/p1g.txt")
p2p = np.genfromtxt("data/p2p.txt")
p2g = np.genfromtxt("data/p2g.txt")

xinds = np.arange(0,19,1)
plt.bar(xinds, p1p, color=martaBlue)
plt.bar(xinds, p1g, bottom=p1p, color=martaRed)
plt.ylim([0,0.5])
plt.show()

plt.clf() 
plt.bar(xinds, p2p, color=martaBlue)
plt.bar(xinds, p2g, bottom=p2p, color=martaRed)
plt.ylim([0,1.65])
plt.show()

#Load in data
# data = np.genfromtxt("data/pixelruns/25_25_0_q6_kap5_D5.txt")
# alpha, S = np.hsplit(data, 2) 

# #Find ac 
# aftertrans = np.argwhere(S > 0.5)
# ac = alpha[aftertrans[0]][0][0]
# print(ac)
# text = "$\\alpha_c = $" + "{:.2f}".format(ac-(2.0/60.0))

# fig, ax = plt.subplots(figsize=(5,4))
# ax.axvline(ac-(2.0/60.0), color="k", lw=3.0, alpha=0.8)
# ax.plot(alpha, S, lw=3.0, color=martaPurple)
# ax.text(1.7, 0.45, text)
# ax.set_xlim([0.5,2.5])
# ax.set_xticks([0.5,2.5])
# ax.set_ylim([0,1])
# ax.set_yticks([0,1])
# ax.set_xlabel("$\\alpha$", labelpad=-10, fontsize=24)
# ax.set_ylabel("$\\mathcal{S}$", rotation=0, fontsize=24)

# plt.tight_layout()
# plt.show() 

# fig = plt.figure()
# ax = fig.add_subplot(111)


# dat = np.genfromtxt("data/swingtimeprofs/q1_10_40_0_svt_alpha125.txt")
# s, t = np.hsplit(dat,2)
# s = s[:,0]
# t = t[:,0]

# xy = np.vstack([t,s])
# z = gaussian_kde(xy)(xy)

# idx = z.argsort()
# t, s, z = t[idx], s[idx], z[idx]

# # fig, ax = plt.subplots()
# ax.scatter(t, s, c=z, s=30, edgecolor='', cmap="magma")
# # ax.set_xlabel("Relative load")
# # ax.set_ylabel("Transient time")
# # ax.set_xlim([0,4])
# # ax.set_ylim([0,20])
# # ax.set_ylabel("Solution distance")
# plt.tight_layout() 
# plt.show()












