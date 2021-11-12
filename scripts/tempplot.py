import numpy as np 
import matplotlib 
import matplotlib.pyplot as plt 
from numpy.random import rand
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

fig = plt.figure()
ax = fig.add_subplot(111)

dat = np.genfromtxt("data/new/q1_30_30_0.txt")
a, s, overs, desy, unb = np.hsplit(dat,5)
a = a[:,0]
s = s[:,0]
overs = overs[:,0]
desy = desy[:,0]
unb = unb[:,0]
# print(np.shape(a))

aftertrans = np.argwhere(s > 0.5)
ac = a[aftertrans[0]][0]
print(ac)

ax.plot(a, s, color=martaBlue, lw=3.0)
ax.axvline(ac, color="k", lw=3.0, alpha=0.7)
# ax.fill_between(a, s*0.0, s, color=martaBlue, alpha=0.2)
ax.plot(a, overs, color=martaGreen, lw=3.0)
ax.plot(a, desy+unb, color=martaRed, lw=3.0)

ax.set_xlim([0.1,2.5])
ax.set_ylim([0,1])
# ax.legend(fontsize=16)
# ax.set_ylabel("$||\\omega-\\omega_0||_2$", labelpad=20)
# ax.set_xlabel("$t$")
plt.tight_layout()
plt.show()












