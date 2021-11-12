import glob, os
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt 
from matplotlib import colors
from ternary.helpers import simplex_iterator
import ternary 
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

#Crawl through the files
search_dir = "data/tri_sm_1/nconthalfhalf/"
files = filter(os.path.isfile, glob.glob(search_dir + "*.txt"))

#Iterate through the files 
fmaxes = []
nstarg, ndtarg, netarg = 50,50,0
cutoff = min([nstarg, ndtarg]) 
nvec = [] 
acvec = []
for i, f in enumerate(files):
	print i, f 

	#Read in file
	dat = np.genfromtxt(f, skip_header=0)

	#Read the config info 
	nst, ndt, net, fmet = np.hsplit(dat[:1], 4)
	ns, nd, ne, fmax = nst[0][0], ndt[0][0], net[0][0], fmet[0][0] 
	n = ns + nd + ne 

	a, e, std, _ = np.hsplit(dat[1:], 4)

	aftertransition = np.argwhere(e > 0.5)
	ac = a[aftertransition[0][0]]
	acvec.append((ac-fmax)/fmax) 
	nvec.append(n)

	# if (ns == nstarg) and (nd == ndtarg) and (ne == netarg):

	# 	#Read the surviving edges vs alpha profile
	# 	a, e, std, power = np.hsplit(dat[1:], 4)

	# 	plt.plot(a, e, ".-", lw=3.0, color=martaRed)
	# 	# plt.xlim([0.2/cutoff,0.5-(cutoff/75.0)])
	# 	plt.show() 

	# 	break 

plt.plot(nvec[:], acvec[:], "o", color=martaBlue)
plt.show() 
