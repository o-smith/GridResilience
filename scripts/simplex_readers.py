import glob, os 
import numpy as np 
import matplotlib 
import matplotlib.pyplot as plt 
import ternary 
from ternary.helpers import simplex_iterator 
from scipy.stats import lognorm 

#Set up plot style
font = {'size'   : 13}
matplotlib.rc('font', **font)
matplotlib.rc('font', serif='Computer Modern Roman')

#Define colours 
martaRed = "#c24c51" 
martaGreen = "#54a666"
martaBlue = "#4c70b0"
martaPurple = "#7f70b0"
martaGold = "#ccb873"


class SimplexReader: 

	def __init__(self, **kwargs):
		self.__dict__.update(kwargs) 
		self.dir = "data/"
		self.critdir = "data/"
		self.n = 50 


	def read_and_store_critcouplings(self, **kwargs):
		self.__dict__.update(kwargs)
		files = filter(os.path.isfile, glob.glob(self.critdir + "*.txt")) 
		print(self.critdir)
		self.crit_matrix = np.zeros((self.n,self.n))
		self.crit_std_matrix = np.zeros((self.n,self.n))
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
			try:
				ind = int(f[41:-4])
			except ValueError: 
				ind = int(f[42:-4])
			ns, nd, ne = self.configs[ind-1][0], self.configs[ind-1][1], self.configs[ind-1][2]
			u, v = int(ns-1), int(nd-1) 
			self.crit_matrix[u,v] = np.mean(dat)
			self.crit_std_matrix[u,v] = np.std(dat)


	def read_and_store_means(self, **kwargs): 
		self.__dict__.update(kwargs) 
		files = filter(os.path.isfile, glob.glob(self.dir + "*.txt")) 
		self.data_matrix = np.zeros((self.n, self.n))
		for i, f in enumerate(files): 
			print(i, f) 
			dat = np.genfromtxt(f, skip_footer=200, dtype=None, invalid_raise=False)
			ns, nd, ne = dat[0], dat[1], dat[2]  
			print(ns, nd)
			dat = np.genfromtxt(f, skip_header=1, dtype=None, invalid_raise=False)
			u, v = int(ns-1), int(nd-1) 
			self.data_matrix[u,v] = np.nanmean(dat)


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


	def plot_means(self, title="simplex", **kwargs): 
		self.__dict__.update(kwargs)
		d = dict() 
		for (i,j,k) in simplex_iterator(self.n-2):
			if self.data_matrix[i,j] == 0.0: 
				d[(j,k)] = np.NaN 
			else: 
				d[(j,k)] = self.data_matrix[i,j]

		figure, tax = ternary.figure(scale=self.n-2)
		figure.set_size_inches(5, 5)
		figure.set_dpi(90)
		tax.boundary(linewidth=0.01)
		tax.gridlines(linewidth=0.000005, multiple=50)
		tax.heatmap(d, style="t", cmap=plt.cm.get_cmap('Reds',10))
		tax.clear_matplotlib_ticks()
		plt.axis('off')
		plt.title(title)
		ternary.plt.show()

	def plot_crits(self, title="simplex", **kwargs): 
		self.__dict__.update(kwargs)
		d = dict() 
		for (i,j,k) in simplex_iterator(self.n-2):
			if self.crit_matrix[i,j] == 0.0: 
				d[(j,k)] = np.NaN 
			else: 
				d[(j,k)] = self.crit_matrix[i,j] 
		figure, tax = ternary.figure(scale=self.n-2)
		figure.set_size_inches(5, 5)
		figure.set_dpi(90)
		tax.boundary(linewidth=0.01)
		tax.gridlines(linewidth=0.000005, multiple=50)
		tax.heatmap(d, style="t", cmap="viridis")
		tax.clear_matplotlib_ticks()
		plt.axis('off')
		plt.title(title)
		ternary.plt.show() 

	def crit_slices(self, **kwargs): 
		self.__dict__.update(kwargs)
		means, stds = [], []
		for (i,j,k) in simplex_iterator(self.n-2): 
			if i + j == 35: 
				means.append(self.crit_matrix[i,j]) 
				stds.append(self.crit_std_matrix[i,j]) 
		return np.array(means), np.array(stds)  



#Read in coupling simplex data
x = SimplexReader()
x.n = 50
x.read_and_store_critcouplings(critdir="data/simplexdata/swing_coupling_q10/50tri/")
means1, stds1 = x.crit_slices()
x.read_and_store_critcouplings(critdir="data/simplexdata/swing_coupling_q4/50tri/")
means2, stds2 = x.crit_slices()
x.read_and_store_critcouplings(critdir="data/simplexdata/swing_coupling_q1/50tri/")
means3, stds3 = x.crit_slices()
x.read_and_store_critcouplings(critdir="data/simplexdata/swing_coupling_q0/50tri/")
means4, stds4 = x.crit_slices()

fig = plt.figure(figsize=(6,4))
ax = fig.add_subplot(111)

#Plot slices through the coupling simplexes
consumers = np.arange(0,36,1)
ax.plot(consumers, means1, lw=3.0, color=martaGold, label="$q=1$")
ax.fill_between(consumers, means1+stds1, means1-stds1, alpha=0.5, color=martaGold)
ax.plot(consumers, means2, lw=3.0, color=martaGreen, label="$q=0.4$")
ax.fill_between(consumers, means2+stds2, means2-stds2, alpha=0.5, color=martaGreen)
ax.plot(consumers, means3, lw=3.0, color=martaBlue, label="$q=0.1$")
ax.fill_between(consumers, means3+stds3, means3-stds3, alpha=0.5, color=martaBlue)
ax.plot(consumers, means4, lw=3.0, color=martaRed, label="$q=0$")
ax.fill_between(consumers, means4+stds4, means4-stds4, alpha=0.5, color=martaRed) 
ax.set_xlim(0,35)
ax.set_title("Mean critical coupling for Watts-Strogatz networks with 50 nodes, 15 of which are passive", fontsize=8)
ax.set_xlabel("Consumers")
ax.set_ylabel("$\\overline{\\kappa}_c$", rotation=0, labelpad=20)

plt.legend()
plt.tight_layout() 
plt.show()

#Plot the whole coupling simplex
print("Plotting critical coupling simplex...")
x.plot_crits("Critical coupling for q=0, n=50")  

#Plot the resilience simplex
x.read_and_store_means(dir="data/simplexdata/swing_rho/50tri/")
print("Plotting rho simplex for n=50 ...")
x.plot_means("$\\overline{\\rho}$ for q=0, n=50") 

#Plot the resilience simplex for n=100 
x.n = 100 
x.read_and_store_means(dir="data/simplexdata/swing_rho/100tri/")
print("Plotting rho simplex for n=100 ...")
x.plot_means("$\\overline{\\rho}$ for q=0, n=100") 






















