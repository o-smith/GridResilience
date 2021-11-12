import glob, os
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt 
from matplotlib import colors
import matplotlib.colors as colors
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


#Crawl through the files
search_dir = "data/qncentre/"
files = filter(os.path.isfile, glob.glob(search_dir + "*.txt"))

#Iterate through the files 
fmaxes = []
datamatrix = np.zeros((40, 40))
x = []  
for i, f in enumerate(files):
	# print(i, f)  

	#Read in file
	dat = np.genfromtxt(f, skip_header=0)
	# print(dat[0])

	#Read the config info  
	i, j, q, n, acc, acm = dat[0], dat[1], dat[2], dat[3], dat[4], dat[5]

	#Find the transition point  
	k = int(i) 
	l = int(j)
	if int(n) == 48:
		x.append(acc - acm)
	# print(type(k))
	datamatrix[k-1,l-1] = acc - acm

# for i in range(51):
# 	for j in range(51):
# 		if datamatrix[i,j] == 0.0:
# 			datamatrix[i,j] = np.NaN
		# if datamatrix[i,j] > 0.0:
		# 	datamatrix[i,j] = 1.0
		# if datamatrix[i,j] < 0.0:
		# 	datamatrix[i,j] = -1.0

# plt.plot(x[::-1]) 
# plt.show() 


#Crawl through the files
search_dir = "data/qncentre2/"
files = filter(os.path.isfile, glob.glob(search_dir + "*.txt"))

#Iterate through the files 
fmaxes = []
datamatrix2 = np.zeros((40, 40))
x2 = []  
for i, f in enumerate(files):
	# print(i, f)  

	#Read in file
	dat = np.genfromtxt(f, skip_header=0)
	# print(dat[0])

	#Read the config info 
	i, j, q, n, acc, acm = dat[0], dat[1], dat[2], dat[3], dat[4], dat[5]

	#Find the transition point 
	k = int(i) 
	l = int(j)
	if int(n) == 120:
		x2.append(acc - acm)
	# print(type(k))
	datamatrix2[k-1, l-1] = acc - acm

#for i in range(51):
# 	for j in range(51):
# 		if datamatrix2[i,j] == 0.0:
# 			datamatrix2[i,j] = np.NaN

y = (datamatrix + datamatrix2)/2.0
# x = (x+x2)/2.0
z = np.transpose(y[:,::-1])
plt.plot(z[-6,:], "o", ms=7, color=martaRed)
plt.ylim([-0.4,0.4])
plt.xlim([0,40])
plt.show()

# for i in range(40):
# 	for j in range(40):
# 		if y[i,j] == 0.0:
# 			y[i,j] = np.NaN
# y = gaussian_filter(y, sigma=0.6)
# cmap=matplotlib.cm.Spectral
# # set the colormap and centre the colorbar
# class MidpointNormalize(colors.Normalize):
# 	"""
# 	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

# 	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
# 	"""
# 	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
# 		self.midpoint = midpoint
# 		colors.Normalize.__init__(self, vmin, vmax, clip)

# 	def __call__(self, value, clip=None):
# 		# I'm ignoring masked values and all kinds of edge cases to make a
# 		# simple example...
# 		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
# 		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

# mid_val = 0.0 
# elev_min = -0.9 #np.min(y)
# elev_max = 0.5
# print(np.min(y))
# plt.imshow(np.transpose(y[:,::-1]), interpolation="bicubic", cmap=cmap, clim=(elev_min, elev_max), norm=MidpointNormalize(midpoint=mid_val,vmin=elev_min, vmax=elev_max),aspect="auto", extent=[0.0,1,40,280])#, vmin=0, vmax=1)
# plt.xlabel("$q$", fontsize=30)
# plt.ylabel("$n$", fontsize=30, rotation=0, labelpad=20)
# plt.colorbar() 
# plt.tight_layout() 
# plt.show() 
