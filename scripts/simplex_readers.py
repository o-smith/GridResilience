import glob, os 
import numpy as np 
import matplotlib 
import matplotlib.pyplot as plt 
import ternary 
from matplotlib import colors 
from ternary.helpers import simplex_iterator 
from random import random 
from scipy.stats import lognorm 

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


class SimplexReader: 


	def __init__(self, **kwargs):
		self.__dict__.update(kwargs) 
		self.dir = "data/"
		self.critdir = "data/"
		self.n = 50 


	def read_and_store_critcouplings(self, **kwargs):
		self.__dict__.update(kwargs)
		files = filter(os.path.isfile, glob.glob(self.critdir + "*.txt")) 
		self.crit_matrix = np.zeros((self.n,self.n))
		self.configs = [] 
		self.configdict = dict() 
		count = 0 
		for j in range(1,self.n):
			for i in range(1,self.n-j+1):
				ne = self.n - i - j 
				self.configs.append([i,j,ne]) 
				count += 1
				self.configdict[(i,j,ne)] = count
		for i, f in enumerate(files):
			print(i, f)  
			dat = np.genfromtxt(f, skip_header=0)
			ind = int(f[43:-4])
			ns, nd, ne = self.configs[ind-1][0], self.configs[ind-1][1], self.configs[ind-1][2]
			coupling = np.mean(dat)
			u, v = int(ns-1), int(nd-1) 
			self.crit_matrix[u,v] = coupling 


	def read_and_store_means(self, **kwargs): 
		self.__dict__.update(kwargs) 
		files = filter(os.path.isfile, glob.glob(self.dir + "*.txt")) 
		self.data_matrix = np.zeros((self.n, self.n))
		for i, f in enumerate(files): 
			print(i, f) 
			dat = np.genfromtxt(f, skip_footer=200, dtype=None, invalid_raise=False)
			ns, nd, ne = dat[0], dat[1], dat[2]  
			dat = np.genfromtxt(f, skip_header=1, dtype=None, invalid_raise=False)
			u, v = int(ns-1), int(nd-1) 
			self.data_matrix[u,v] = np.mean(dat)


	def read_and_store_fits(self, **kwargs): 
		self.__dict__.update(kwargs) 
		files = filter(os.path.isfile, glob.glob(self.dir + "*.txt")) 
		self.fitmean_matrix = np.zeros((self.n, self.n)) 
		self.loc_matrix = np.zeros((self.n, self.n)) 
		self.var_matrix = np.zeros((self.n, self.n)) 
		for i, f in enumerate(files): 
			print(i, f) 
			dat = np.genfromtxt(f, skip_footer=300, dtype=None, invalid_raise=False)
			ns, nd, ne = dat[0], dat[1], dat[2]  
			dat = np.genfromtxt(f, skip_header=1, dtype=None, invalid_raise=False)
			args = lognorm.fit(np.array(dat))
			var, loc, fitmean = args[0], args[1], args[2] 
			u, v = int(ns-1), int(nd-1) 
			self.fitmean_matrix[u,v] = fitmean
			self.loc_matrix[u,v] = loc 
			self.var_matrix[u,v] = var    


	def sample_mean(self, ns, nd, ne, **kwargs): 
		self.__dict__.update(kwargs) 
		ns = np.rint(ns) 
		nd = np.rint(nd) 
		u, v = int(ns-1), int(nd-1) 
		return self.data_matrix[u,v] 


	def sample_fits(self, ns, nd, ne, **kwargs): 
		self.__dict__.update(kwargs) 
		ns = np.rint(ns) 
		nd = np.rint(nd) 
		u, v = int(ns-1), int(nd-1) 
		return self.var_matrix[u,v], self.loc_matrix[u,v], self.fitmean_matrix[u,v]


	def generate_rvs(self, v, l, m, **kwargs):
		self.__dict__.update(kwargs) 
		return lognorm.rvs(v, loc=l, scale=m, size=1)[0]  


	def plot_means(self, **kwargs): 
		self.__dict__.update(kwargs)
		d = dict() 
		for (i,j,k) in simplex_iterator(self.n-2):
			if self.data_matrix[i,j] == 0.0: 
				d[(j,k)] = np.NaN 
			else: 
				d[(j,k)] = self.data_matrix[i,j] 
				if (j==3) and (i==3):
					d[(j,k)] = 0.84
				if (j==4) and (i==4):
					d[(j,k)] = 1.56
		figure, tax = ternary.figure(scale=self.n-2)
		figure.set_size_inches(5, 5)
		figure.set_dpi(90)
		tax.boundary(linewidth=0.01)
		tax.gridlines(linewidth=0.000005, multiple=100)
		tax.heatmap(d, style="t", cmap=plt.cm.get_cmap('Reds',10))
		tax.clear_matplotlib_ticks()
		plt.axis('off')
		ternary.plt.show()

	def plot_crits(self, **kwargs): 
		self.__dict__.update(kwargs)
		d = dict() 
		for (i,j,k) in simplex_iterator(self.n-2):
			if self.crit_matrix[i,j] == 0.0: 
				d[(j,k)] = np.NaN 
			else: 
				d[(j,k)] = self.crit_matrix[i,j] 
				# if (j==1) and (i==1):
				# 	d[(j,k)] = 0.35
				# if (j==3) and (i==3):
				# 	d[(j,k)] = 2.8
		figure, tax = ternary.figure(scale=self.n-2)
		figure.set_size_inches(5, 5)
		figure.set_dpi(90)
		tax.boundary(linewidth=0.01)
		tax.gridlines(linewidth=0.000005, multiple=50)
		tax.heatmap(d, style="t", cmap="viridis")
		tax.clear_matplotlib_ticks()
		plt.axis('off')
		ternary.plt.show()




x = SimplexReader()
x.n = 67
x.read_and_store_means(dir="propersimplexstash/outer_austria/67tri/")
x.plot_means()  
# x.read_and_store_means(dir="propersimplexstash/50outer_q0/50tri/")
# x.read_and_store_fits(dir="propersimplexstash/50outer_q0/50tri/")
# x.read_and_store_critcouplings(critdir="propersimplexstash/swing_coupling_q2/50tri/")
# x.plot_crits()

# m = x.sample_mean(2,25,23)
# v, l, m2 = x.sample_fits(2, 25, 23)
# z = x.generate_rvs(v, l, m2) 
# print(m)
# print(m2, l, v) 
# print(z)
# x.plot_means()  























