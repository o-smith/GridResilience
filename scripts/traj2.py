import glob, os
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt 
from matplotlib import colors
from ternary.helpers import simplex_iterator
import ternary 
from scipy.optimize import curve_fit
from scipy.ndimage.filters import gaussian_filter
from random import random

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

#Read in simplex data 
#Iterate through the files
#Crawl through the files
dim = 50
search_dir = "data/50tri/"
files = filter(os.path.isfile, glob.glob(search_dir + "*.txt")) 
fmaxes = []
datamatrix = np.zeros((dim, dim)) 
for i, f in enumerate(files):
	print(i, f) 

	#Read in file
	dat = np.genfromtxt(f, skip_header=0)

	#Read the config info 
	nst, ndt, net, fmet = np.hsplit(dat[:1], 4)
	ns, nd, ne, fmax = nst[0][0], ndt[0][0], net[0][0], fmet[0][0] 

	#Read the surviving edges vs alpha profile
	ac, fmax, iters, iters = np.hsplit(dat[1:], 4)

	u, v = int(ns-1), int(nd-1)
	#fmaxes.append(fmax)
	datamatrix[u,v] = ac


#Read in trajectory data
trace = np.genfromtxt("data/adayinspring4.txt")
remapped_configs = [] 
i = 0
for config in trace:
	ns = config[0]
	nd = config[1]
	ne = config[2] 
	x = (nd, ne, ns)
	remapped_configs.append(x)

resilience = [] 
for element in remapped_configs: 
	nd, ne, ns = element
	u, v = int(ns-1), int(nd-1)
	resi = datamatrix[u,v]
	resilience.append(resi)
ax.plot(resilience)

# #Read in trajectory data
# trace = np.genfromtxt("data/adayinspring2.txt")
# remapped_configs = [] 
# i = 0
# for config in trace:
# 	ns = config[0]
# 	nd = config[1]
# 	ne = config[2] 
# 	x = (nd, ne, ns)
# 	remapped_configs.append(x)

# resilience = [] 
# for element in remapped_configs: 
# 	nd, ne, ns = element
# 	u, v = int(ns-1), int(nd-1)
# 	resi = datamatrix[u,v]
# 	resilience.append(resi)
# ax.plot(resilience)

# #Read in trajectory data
# trace = np.genfromtxt("data/adayinspring3.txt")
# remapped_configs = [] 
# i = 0
# for config in trace:
# 	ns = config[0]
# 	nd = config[1]
# 	ne = config[2] 
# 	x = (nd, ne, ns)
# 	remapped_configs.append(x)

# resilience = [] 
# for element in remapped_configs: 
# 	nd, ne, ns = element
# 	u, v = int(ns-1), int(nd-1)
# 	resi = datamatrix[u,v]
# 	resilience.append(resi)
# ax.plot(resilience)

# #Read in trajectory data
# trace = np.genfromtxt("data/adayinspring.txt")
# remapped_configs = [] 
# i = 0
# for config in trace:
# 	ns = config[0]
# 	nd = config[1]
# 	ne = config[2] 
# 	x = (nd, ne, ns)
# 	remapped_configs.append(x)

# resilience = [] 
# for element in remapped_configs: 
# 	nd, ne, ns = element
# 	u, v = int(ns-1), int(nd-1)
# 	resi = datamatrix[u,v]
# 	resilience.append(resi)
# ax.plot(resilience)
# ax.set_xticks([])
# ax.set_ylabel("$\\rho$", rotation=0)

plt.tight_layout() 
plt.show() 










