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

dim = 40

#Set up ternary plot
scale = dim 
figure, tax = ternary.figure(scale=scale)
figure.set_size_inches(6, 5.5)
figure.set_dpi(90)

# Draw Boundary and grid lines
tax.boundary(linewidth=1)

#Set Axis labels 
fontsize = 26
tax.left_axis_label("$\\leftarrow $ Source nodes" , fontsize=fontsize)
tax.right_axis_label("$\\leftarrow$ Passive nodes", fontsize=fontsize)
tax.bottom_axis_label("Sink nodes $\\rightarrow $", fontsize=fontsize)


p1 = (22, 8, 10)
p2 = (2, 22, 16)

# trace = np.genfromtxt("data/adayinspring4.txt")

# remapped_configs = [] 
# i = 0
# for config in trace:
# 	ns = config[0]
# 	nd = config[1]
# 	ne = config[2] 
# 	x = (nd, ne, ns)
# 	remapped_configs.append(x)

	# figure, tax = ternary.figure(scale=scale)
	# figure.set_size_inches(6, 5.5)
	# figure.set_dpi(90)

	# # Draw Boundary and grid lines
	# tax.boundary(linewidth=1)

	# #Set Axis labels 
	# fontsize = 26
	# tax.left_axis_label("$\\leftarrow $ Source nodes" , fontsize=fontsize)
	# tax.right_axis_label("$\\leftarrow$ Passive nodes", fontsize=fontsize)
	# tax.bottom_axis_label("Sink nodes $\\rightarrow $", fontsize=fontsize)
	# tax.gridlines(linewidth=0.5, multiple=2)
	# tax.plot(remapped_configs, color=martaRed)
	# tax.scatter([x], color=martaRed)
	# # tax.ticks(axis='lbr', multiple=10, linewidth=2, offset=0.025)
	# tax.get_axes().axis('off')
	# tax.clear_matplotlib_ticks()
	# fname = "frame%i.png" %i
	# tax.savefig(fname)
	# plt.clf()
	# plt.cla()
	# i += 1
# trace = np.genfromtxt("data/adayinspring4.txt")

# remapped_configs = [] 
# i = 0
# for config in trace:
# 	ns = config[0]
# 	nd = config[1]
# 	ne = config[2] 
# 	x = (nd, ne, ns)
# 	remapped_configs.append(x)

# tax.gridlines(linewidth=0.5, multiple=2)
# tax.plot(remapped_configs, color="k")

# trace = np.genfromtxt("data/adayinspring5.txt")

# remapped_configs = [] 
# i = 0
# for config in trace:
# 	ns = config[0]
# 	nd = config[1]
# 	ne = config[2] 
# 	x = (nd, ne, ns)
# 	remapped_configs.append(x)

# tax.gridlines(linewidth=0.5, multiple=2)
# tax.plot(remapped_configs, color="c")

# trace = np.genfromtxt("data/adayinspring3.txt")

# remapped_configs = [] 
# i = 0
# for config in trace:
# 	ns = config[0]
# 	nd = config[1]
# 	ne = config[2] 
# 	x = (nd, ne, ns)
# 	remapped_configs.append(x)

# tax.gridlines(linewidth=0.5, multiple=2)
# tax.plot(remapped_configs, color="k")

# trace = np.genfromtxt("data/adayinspring2.txt")

# remapped_configs = [] 
# i = 0
# for config in trace:
# 	ns = config[0]
# 	nd = config[1]
# 	ne = config[2] 
# 	x = (nd, ne, ns)
# 	remapped_configs.append(x)

# tax.gridlines(linewidth=0.5, multiple=2)
# tax.plot(remapped_configs, color="r")

# trace = np.genfromtxt("data/adayinspring.txt")

# remapped_configs = [] 
# i = 0
# for config in trace:
# 	ns = config[0]
# 	nd = config[1]
# 	ne = config[2] 
# 	x = (nd, ne, ns)
# 	remapped_configs.append(x)

# tax.gridlines(linewidth=0.5, multiple=2)
# tax.plot(remapped_configs, color="g")

trace = np.genfromtxt("data/adayinspring6.txt")

remapped_configs = [] 
i = 0
for config in trace:
	ns = config[0]
	nd = config[1]
	ne = config[2] 
	x = (nd, ne, ns)
	remapped_configs.append(x)

tax.gridlines(linewidth=0.5, multiple=2)
tax.plot(remapped_configs, linewidth=3.0, color=martaBlue)
# tax.scatter(remapped_configs)
# tax.ticks(axis='lbr', multiple=10, linewidth=2, offset=0.025)
tax.get_axes().axis('off')
tax.clear_matplotlib_ticks()
tax.show()


