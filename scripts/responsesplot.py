import numpy as np 
import matplotlib 
import matplotlib.pyplot as plt 

#Set up plot style
font = {'size'   : 12}
matplotlib.rc('font', **font)
matplotlib.rc('font', serif='Computer Modern Roman')

#Define colours 
martaRed = "#c24c51" 
martaBlue = "#4c70b0"
D4colour = '#97015E' #Metropolitan line magenta


fig = plt.figure()
ax2 = fig.add_subplot(222) 

Ds = np.arange(3,20,1)
alphas = np.linspace(0.3,1.0,len(Ds))
for i, D in enumerate(Ds):
	print(alphas[i])
	fname = "data/durations/q1_10_40_0_t_gamma" + str(int(D)) + ".txt"
	dat = np.genfromtxt(fname)
	a, t = np.hsplit(dat, 2) 
	ax2.plot(a, t, alpha=alphas[i], color=martaBlue)

kappas = np.arange(1,12,1)
alphas = np.linspace(0.2,1.0,len(kappas))
ax3 = fig.add_subplot(223) 
for i, k in enumerate(kappas):
	print(alphas[i])
	fname = "data/durations/q1_10_40_0_t_kappa" + str(int(k)) + ".txt"
	dat = np.genfromtxt(fname)
	a, t = np.hsplit(dat, 2) 
	ax3.plot(a, t, alpha=alphas[i], color="k")

ns = np.arange(45,105,5)
ax1 = fig.add_subplot(221) 
alphas = np.linspace(0.3,1.0,len(ns))
for i, n in enumerate(ns):
	print(alphas[i])
	fname = "data/durations/q1_t_" + str(int(n)) + ".txt"
	dat = np.genfromtxt(fname)
	a, t = np.hsplit(dat, 2) 
	ax1.plot(a, t, alpha=alphas[i], color=martaRed)

qs = np.arange(0,11,1)
ax4 = fig.add_subplot(224) 
alphas = np.linspace(0.2,1.0,len(qs))
for i, q in enumerate(qs):
	print(alphas[i])
	fname = "data/durations/q1_10_40_0_t_q" + str(int(q)) + ".txt"
	dat = np.genfromtxt(fname)
	a, t = np.hsplit(dat, 2) 
	ax4.plot(a, t, alpha=alphas[i], color=D4colour)


ax1.set_xlim([0.5,2.5])
ax2.set_xlim([0.5,2.5])
ax3.set_xlim([0.5,2.5])
ax4.set_xlim([0.5,2.5])
ax1.set_title("Increasing $n$")
ax2.set_title("Increasing $\\gamma$")
ax3.set_title("Increasing $\\kappa$")
ax4.set_title("Increasing $q$")
ax1.set_xlabel("$\\alpha/\\alpha_{\\ast}$")
ax1.set_ylabel("$\\overline{T}$", rotation=0, labelpad=5)
ax2.set_xlabel("$\\alpha/\\alpha_{\\ast}$")
ax2.set_ylabel("$\\overline{T}$", rotation=0, labelpad=5)
ax3.set_xlabel("$\\alpha/\\alpha_{\\ast}$")
ax3.set_ylabel("$\\overline{T}$", rotation=0, labelpad=5)
ax4.set_xlabel("$\\alpha/\\alpha_{\\ast}$")
ax4.set_ylabel("$\\overline{T}$", rotation=0, labelpad=5)


plt.tight_layout()
plt.show()

