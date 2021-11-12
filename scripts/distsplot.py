import glob, os
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt 
from matplotlib import colors
from mpl_toolkits.mplot3d import Axes3D
# from ternary.helpers import simplex_iterator
# import ternary 
from scipy.optimize import curve_fit
from scipy.ndimage.filters import gaussian_filter
from scipy.stats import lognorm, gamma, kstest, expon  

#Set up plot style
font = {'size'   : 18}
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


fig, ax = plt.subplots(figsize=(5,4))
da = np.genfromtxt("data/new/hist_q1_30_30_0.txt")
s, dat = np.hsplit(da,2)
# print(np.mean(dataraw))
# dat = list(filter(lambda x: (x<2.5 and np.isnan(x)==False), dataraw))
H, bins = np.histogram(np.array(dat), density=True, bins=30)
xs = (bins[:-1] + bins[1:])/2
# ax.fill_between(xs/np.mean(dat), H*np.mean(dat)) #, zs=0.0, zdir='y')
ax.bar(xs/np.mean(dat), H*np.mean(dat), width=max(dat)/30.0,facecolor=martaRed, alpha=0.7)
args = lognorm.fit(np.array(dat)) 
# print(args)
ks, pval = kstest(dat, 'lognorm', args=args)
print("kstest: ", ks, "pval: ", pval)
pdf_exp = lognorm.pdf(bins, *args)
ax.plot(bins/np.mean(dat), pdf_exp*np.mean(dat), ls='-', color="k", lw=3.2)
# ax.set_xlabel("$\\kappa_c/\\bar{\\kappa}_c$")
# ax.set_ylabel("$\\bar{\\kappa}_c P(\\kappa_c)$", rotation=0, labelpad=20)
# ax.set_xlim([0.35,1.8])
# ax.set_xticks([0.35,1.8])
ax.set_yticks([0,3])
ax.set_ylim([0,3])
plt.tight_layout()
plt.show()


# fig = plt.figure(figsize=(6,3.3))
# ax = fig.add_subplot(111, projection='3d')


# # ax1 = fig.add_subplot(121)
# # ax2 = fig.add_subplot(122)

# fonstsz = 13

# dat = np.genfromtxt("data/achists/sm_30_30_q1_new.txt")
# H, bins = np.histogram(np.array(dat), density=True, bins=40)
# xs = (bins[:-1] + bins[1:])/2
# # ax.fill_between(xs/np.mean(dat), H*np.mean(dat)) #, zs=0.0, zdir='y')
# ax.bar(xs/np.mean(dat), H*np.mean(dat), zs=0.0, zdir='y', width=max(dat)/70.0,facecolor=martaRed, alpha=0.7)
# args = lognorm.fit(np.array(dat), loc=0, scale=1) 
# pdf_exp = lognorm.pdf(bins, *args)
# ax.plot(bins/np.mean(dat), pdf_exp*np.mean(dat), zs=0.0, zdir='y', ls='-', color="k", lw=3.2)
# # ax.plot(bins/np.mean(dat), pdf_exp*np.mean(dat), zs=0.0, zdir='y', ls='-', color=martaRed, lw=3.5) 


# dat = np.genfromtxt("data/achists/sm_30_30_q6_new.txt")
# H, bins = np.histogram(np.array(dat), density=True, bins=40)
# xs = (bins[:-1] + bins[1:])/2
# ax.bar(xs/np.mean(dat), H*np.mean(dat), zs=0.3, zdir='y', width=max(dat)/70.0,facecolor=martaBlue, alpha=0.7)
# args = lognorm.fit(np.array(dat), loc=0, scale=1) 
# pdf_exp = lognorm.pdf(bins, *args)
# ax.plot(bins/np.mean(dat), pdf_exp*np.mean(dat), zs=0.3, zdir='y', ls='-', color="k", lw=3.2)
# # ax.plot(bins/np.mean(dat), pdf_exp*np.mean(dat), zs=0.2, zdir='y', ls='-', color=martaBlue, lw=3.5) 

# dat = np.genfromtxt("data/achists/sf_30_30_k4_new.txt")

# #First plot Nash flow distributions
# H, bins2 = np.histogram(np.array(dat), density=True, bins=40)
# ax.bar(bins2[:-1]/np.mean(dat), H*np.mean(dat), zs=0.8, zdir='y', width=max(dat)/70.0,facecolor=martaPurple, alpha=0.7)


# dat = np.genfromtxt("data/achists/sf_30_30_k3_new.txt")

# #First plot Nash flow distributions
# H, bins2 = np.histogram(np.array(dat), density=True, bins=40)
# ax.bar(bins2[:-1]/np.mean(dat), H*np.mean(dat), zs=1.1, zdir='y', width=max(dat)/70.0,facecolor=martaGreen, alpha=0.7)
# # ax.set_ylabel("$\\overline{ \\rho} P(\\rho)$", rotation=0, fontsize=20, labelpad=20)
# # ax.set_xlabel("$\\rho/\\overline{\\rho}$", fontsize=20, labelpad=10)
# # ax.legend(loc=2, prop={'size': 10})
# # ax.text(-0.2, 1.2, "(b)",transform=ax2.transAxes,
#        # va='top', ha='right')

# # #First plot Nash flow distributions
# # H, bins2 = np.histogram(np.array(dat), density=True, bins=50)
# # ax1.bar(bins2[:-1]/np.mean(dat), H*np.mean(dat), width=max(dat)/45.0,facecolor=martaRed, alpha=0.5, label="$q=0.1$")
# # args = lognorm.fit(np.array(dat), loc=0, scale=1) 
# # pdf_exp = lognorm.pdf(bins2, *args)
# # print(args)
# # ax1.plot(bins2/np.mean(dat), pdf_exp*np.mean(dat), '-', color="k", lw=4.2)
# # ax1.plot(bins2/np.mean(dat), pdf_exp*np.mean(dat), '-', color=martaRed, lw=3.5) 

# # ks, pval = kstest(dat, 'lognorm', args=args)
# # print(ks, pval)

# # dat = np.genfromtxt("data/achists/sm_30_30_q6_new.txt")

# # #First plot Nash flow distributions
# # H, bins2 = np.histogram(np.array(dat), density=True, bins=50)
# # ax1.bar(bins2[:-1]/np.mean(dat), H*np.mean(dat), width=max(dat)/45.0,facecolor=martaBlue, alpha=0.6, label="$q=0.6$")
# # args = lognorm.fit(np.array(dat), loc=0, scale=1) 
# # pdf_exp = lognorm.pdf(bins2, *args) 
# # print(args)
# # ax1.plot(bins2/np.mean(dat), pdf_exp*np.mean(dat), '-', color="k", lw=4.2) 
# # ax1.plot(bins2/np.mean(dat), pdf_exp*np.mean(dat), '-', color=martaBlue, lw=3.5) 
# # ax1.set_xlim([0,2.5])
# # ax1.set_xticks([0,2.5])
# # ax1.set_ylim([0,2])
# # ax1.set_yticks([0,2])
# # ax1.set_ylabel("$\\overline{ \\rho} P(\\rho)$", rotation=0, fontsize=20, labelpad=20)
# # ax1.set_xlabel("$\\rho/\\overline{\\rho}$", fontsize=20, labelpad=10)
# # ax1.legend(prop={'size': 10})
# # ax1.text(-0.2, 1.2, "(a)",transform=ax1.transAxes,
# #        va='top', ha='right')
# # ks, pval = kstest(dat, 'lognorm', args=args)
# # print(ks, pval)

# # dat = np.genfromtxt("data/achists/sf_30_30_k3_new.txt")

# # #First plot Nash flow distributions
# # H, bins2 = np.histogram(np.array(dat), density=True, bins=50)
# # ax2.bar(bins2[:-1]/np.mean(dat), H*np.mean(dat), width=max(dat)/45.0,facecolor=martaRed, alpha=0.5, label="$m_0=3$")


# # dat = np.genfromtxt("data/achists/sf_30_30_k4_new.txt")

# # #First plot Nash flow distributions
# # H, bins2 = np.histogram(np.array(dat), density=True, bins=50)
# # ax2.bar(bins2[:-1]/np.mean(dat), H*np.mean(dat), width=max(dat)/45.0,facecolor=martaBlue, alpha=0.5, label="$m_0=4$")
# # ax2.set_xlim([0,1.4])
# # ax2.set_xticks([0,1.4])
# # ax2.set_ylim([0,6])
# # ax2.set_yticks([0,6])
# # ax2.set_ylabel("$\\overline{ \\rho} P(\\rho)$", rotation=0, fontsize=20, labelpad=20)
# # ax2.set_xlabel("$\\rho/\\overline{\\rho}$", fontsize=20, labelpad=10)
# # ax2.legend(loc=2, prop={'size': 10})
# # ax2.text(-0.2, 1.2, "(b)",transform=ax2.transAxes,
# #        va='top', ha='right')

# # plt.tight_layout() 
# # plt.subplots_adjust(left=0.2, wspace=0.5, hspace=0.6)
# # ax.set_ylim([0,0.75])
# # ax.set_xlim([0.4,2])
# ax.set_yticks([])
# ax.set_zlim([0,3.5])
# # ax.grid(False
# plt.show() 









