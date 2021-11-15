import numpy as np 
import matplotlib 
import matplotlib.pyplot as plt 

#Set up plot style
font = {'size'   : 14}
matplotlib.rc('font', **font)
matplotlib.rc('font', serif='Computer Modern Roman')

#Define colours 
martaRed = "#c24c51" 
martaGreen = "#54a666"
martaBlue = "#4c70b0"
martaPurple = "#7f70b0"
martaGold = "#ccb873"

fig = plt.figure(figsize=(12,4))

#Plot cascade results as a function of alpha for (15,45,0) case
ax1 = fig.add_subplot(131)
dat = np.genfromtxt("data/cascadefailurebreakdown1.txt")
a, s, ov, de = np.hsplit(dat, 4) 
ax1.plot(a, s, lw=2, color=martaBlue, label="Survivors")
ax1.plot(a, ov, lw=2, color=martaGreen, label="Overloads")
ax1.plot(a, de, lw=2, color=martaRed, label="De-syncs") 
ax1.set_title("$(n_+,n_-,n_p)=(15,45,0)$") 
ax1.set_ylabel("Proportion", labelpad=10)
ax1.set_xlabel("$\\alpha/\\alpha_{\\ast}$")
ax1.set_ylim([0,1])
ax1.set_xlim([0.1,2.5]) 
ax1.legend()

#Plot cascade results as a function of alpha for (30,30,0) case
ax2 = fig.add_subplot(132)
dat = np.genfromtxt("data/cascadefailurebreakdown2.txt")
a, s, ov, de = np.hsplit(dat, 4) 
ax2.plot(a, s, lw=2, color=martaBlue, label="Survivors")
ax2.plot(a, ov, lw=2, color=martaGreen, label="Overloads")
ax2.plot(a, de, lw=2, color=martaRed, label="De-syncs") 
ax2.set_title("$(n_+,n_-,n_p)=(30,30,0)$") 
ax2.set_ylabel("Proportion", labelpad=10)
ax2.set_xlabel("$\\alpha/\\alpha_{\\ast}$")
ax2.set_ylim([0,1])
ax2.set_xlim([0.1,2.5]) 
ax2.legend() 

#Now plot mean cascade durations 
ax3 = fig.add_subplot(133)
dat = np.genfromtxt("data/cascadedurations1.txt")
a, t = np.hsplit(dat, 2) 
ax3.plot(a, t, lw=2, color="k", alpha=0.7, label="$(n_+,n_-,n_p)=(15,45,0)$")
dat = np.genfromtxt("data/cascadedurations2.txt")
a, t = np.hsplit(dat, 2) 
ax3.plot(a, t, lw=2, color="k", label="$(n_+,n_-,n_p)=(30,30,0)$")
ax3.set_ylim([0,30])
ax3.set_xlim([0.1,2.5])
ax3.set_title("Mean durations of cascades")
ax3.set_ylabel("Mean time (s)", labelpad=10)
ax3.set_xlabel("$\\alpha/\\alpha_{\\ast}$")
ax3.legend() 

plt.tight_layout()
plt.show() 










