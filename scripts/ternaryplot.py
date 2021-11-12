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

dim = 100
n = dim

#Set up ternary plot
scale = dim-2 
figure, tax = ternary.figure(scale=scale)
figure.set_size_inches(5, 5)
figure.set_dpi(90)

# Draw Boundary and grid lines
tax.boundary(linewidth=0.01)

#Set Axis labels 
fontsize = 26
tax.left_axis_label("$\\leftarrow $ Source nodes" , fontsize=fontsize)
tax.right_axis_label("$\\leftarrow$ Passive nodes", fontsize=fontsize)
tax.bottom_axis_label("Sink nodes $\\rightarrow $", fontsize=fontsize)

#Crawl through the files
search_dir = "data/swingq10_I1_D1_kap5/100tri/"
files = filter(os.path.isfile, glob.glob(search_dir + "*.txt"))

configdict = dict() 

# configs = [] 
# n = 67
# count = 0 
# for j in range(1,n):
# 	for i in range(1,n-j+1):
# 		ne = n - i - j 
# 		configs.append([i,j,ne]) 
# 		count += 1
# 		configdict[(i,j,ne)] = count 
		

#Iterate through the files 
fmaxes = []
datamatrix = np.zeros((dim, dim)) 
stdmatrix = np.zeros((dim, dim))
for i, f in enumerate(files):
	print(i, f) 

	#Read in file
	dat = np.genfromtxt(f, skip_header=0)

	#Get index 

	# ind = int(f[23:-4])
	# print(ind) 

	#Get config 
	# ns, nd, ne = configs[ind-1][0], configs[ind-1][1], configs[ind-1][2]

	#Read the config info 
	nst, ndt, net, fmet = np.hsplit(dat[:1], 4)
	ns, nd, ne, fmax = nst[0][0], ndt[0][0], net[0][0], fmet[0][0] 

	#Read the surviving edges vs alpha profile
	ac, fmax, iters, iters = np.hsplit(dat[1:], 4)
	print(ac)
	# coupling = np.mean(dat)

	#Find the transition point  
	#aftertransition = np.argwhere(e > 0.5)
	#ac = a[aftertransition[0][0]]

	u, v = int(ns-1), int(nd-1)
	#fmaxes.append(fmax)
	datamatrix[u,v] = ac
	# datamatrix[u,v] = coupling
	stdmatrix[u,v] = np.std(dat) 

#Copy the matrix into a dictionary
missing_indices = [] 
midslice = []
d = dict() 
for (i,j,k) in simplex_iterator(dim-2):
	if datamatrix[i,j] == 0.0:
		d[(j,k)] = np.NaN
		# print(i,j) 
		# ns = i + 1
		# nd = j + 1 
		# ne = n - ns - nd 
		# # offending_count = configdict[(ns, nd, ne)]
		# missing_indices.append(offending_count)
	else:
		d[(j,k)] = datamatrix[i,j] 
		# ns = i + 1
		# nd = j + 1 
		# ne = n - ns - nd 
		# if ne ==30:
		# 	d[(j,k)] = 0.0
		# 	midslice.append((ns,nd,datamatrix[i,j],stdmatrix[i,j]))
		if (j==10) and (i==10):
			d[(j,k)] = 0.35
		if (j==10) and (i==11):
			d[(j,k)] = 2.8


# # print(missing_indices) 
# sortedslice = list(sorted(midslice, key=lambda tup: tup[1]))
# print(sortedslice)
# outvec = []
# outvec2 = []
# for i in sortedslice:
# 	outvec.append(i[2])
# 	outvec2.append(i[3])
# outvec = np.array(outvec)
# outvec2 = np.array(outvec2) 
# v = zip(outvec, outvec2)
# np.savetxt("q10_slice.txt", v)

# p1 = (0, 148, 0)
# p2 = (75, 0, 75)
# tax.line(p1, p2, linewidth=8., marker='None', alpha=0.8, color=martaBlue, linestyle="--")

# tax.horizontal_line(60, linewidth=8., color=martaGreen, alpha=0.9, linestyle="--")
# tax.left_parallel_line(0, linewidth=8., color='k', alpha=0.6, linestyle='--')

#Plot as a heat map
tax.gridlines(linewidth=0.000005, multiple=100)
tax.heatmap(d, style="t", cmap=plt.cm.get_cmap('Reds',7))
tax.clear_matplotlib_ticks()
plt.axis('off')
# plt.text(46,127,"(iii)", fontsize=28)
# plt.text(40,40,"(ii)", fontsize=28)
# plt.text(80,6,"(i)", fontsize=28)
ternary.plt.show()


# #Set half of the matix to nan  
# vec = [] 
# vec2 = []
# width = np.shape(datamatrix)[0]
# for i in range(width):
# 	for j in range(width):
# 		if datamatrix[i,j] == 0.0:
# 			datamatrix[i,j] = np.NaN 
# 		if i == 1:
# 			if datamatrix[i,j] == np.nan:
# 				continue
# 			else:
# 				vec.append(datamatrix[i,j])
# 				vec2.append(j+1.0)

# # Function to fit to
# func = lambda x, a, b: a/np.sqrt(x) + b
# # func = lambda x, a, b, c: (a*x + b)**2  + c

# print(vec)
# # print(vec2[::-1])
# vec = np.array(vec)[:58]
# vec2 = np.array(vec2)[:58]
# # vec2 = vec2[::-1]

# p0 = np.zeros(2) + 1.0 
# popt, pconv = curve_fit(func, vec2, vec)



# print(popt)
# # perr = np.sqrt(np.diag(pconv))
# # print(perr)
# plt.plot(vec2, vec, 'o', color=martaBlue, alpha=0.6, ms=9.0)
# plt.plot(vec2, func(vec2, *popt), '-', color=martaGreen, lw=5.0, label="$\\frac{0.19}{n_d^2}+0.17$")
# plt.xticks([1,30], fontsize=24)
# # plt.xlim([1,75])
# plt.yticks([0.15,0.45], fontsize=24)
# plt.xlabel("$n_d$", fontsize=34, labelpad=-20)
# plt.ylabel("$\\left<{a_c}\\right>$", fontsize=34, labelpad=-20, rotation=0)
# plt.legend()
# plt.tight_layout() 
# plt.show() 
















