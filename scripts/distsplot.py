import glob, os
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt 
from scipy.stats import lognorm, kstest 

#Set up plot style
font = {'size'   : 18}
matplotlib.rc('font', **font)
matplotlib.rc('font', serif='Computer Modern Roman')

#Define colours 
martaRed = "#c24c51" 


fig = plt.figure(figsize=(11,4))
ax = fig.add_subplot(121)
da = np.genfromtxt("data/hist_q0_50_50_0.txt")
s, dat = np.hsplit(da,2)
H, bins = np.histogram(np.array(dat), density=True, bins=30)
xs = (bins[:-1] + bins[1:])/2
ax.bar(xs/np.mean(dat), H*np.mean(dat), width=max(dat)/30.0,facecolor=martaRed, alpha=0.7)
args = lognorm.fit(np.array(dat)) 
ks, pval = kstest(dat, 'lognorm', args=args)
print("Fit for q=0 case:")
print("kstest: ", ks, "pval: ", pval)
pdf_exp = lognorm.pdf(bins, *args)
ax.plot(bins/np.mean(dat), pdf_exp*np.mean(dat), ls='-', color="k", lw=3.2)
ax.set_title("$\\rho$ dist ofr $q=0$")
ax.set_xlabel("$\\rho/\\bar{\\rho}$")
ax.set_ylabel("$\\bar{\\rho} P(\\rho)$", rotation=0, labelpad=20)
ax.set_xlim([0.35,1.9])
ax.set_xticks([0.35,1.9])
ax.set_yticks([0,3])
ax.set_ylim([0,3])


ax = fig.add_subplot(122)
da = np.genfromtxt("data/hist_q1_50_50_0.txt")
s, dat = np.hsplit(da,2)
H, bins = np.histogram(np.array(dat), density=True, bins=30)
xs = (bins[:-1] + bins[1:])/2
ax.bar(xs/np.mean(dat), H*np.mean(dat), width=max(dat)/30.0,facecolor=martaRed, alpha=0.7)
args = lognorm.fit(np.array(dat)) 
ks, pval = kstest(dat, 'lognorm', args=args)
print("Fit for q=0.1 case:")
print("kstest: ", ks, "pval: ", pval)
pdf_exp = lognorm.pdf(bins, *args)
ax.plot(bins/np.mean(dat), pdf_exp*np.mean(dat), ls='-', color="k", lw=3.2)
ax.set_title("$\\rho$ dist ofr $q=0.1$")
ax.set_xlabel("$\\rho/\\bar{\\rho}$")
ax.set_ylabel("$\\bar{\\rho} P(\\rho)$", rotation=0, labelpad=20)
ax.set_xlim([0.35,1.9])
ax.set_xticks([0.35,1.9])
ax.set_yticks([0,3])
ax.set_ylim([0,3])


plt.tight_layout()
plt.show()









