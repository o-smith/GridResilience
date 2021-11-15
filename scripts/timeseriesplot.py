import matplotlib
from powerreader import * 

#Set up plot style
font = {'size'   : 12}
matplotlib.rc('font', **font)
matplotlib.rc('font', serif='Computer Modern Roman')

#Define colours 
martaRed = "#c24c51" 
martaGreen = "#54a666"
martaBlue = "#4c70b0"
martaPurple = "#7f70b0"
martaGold = "#ccb873"

dates, stds, means, ns = wholeset_meanyears_PV() 

fig = plt.figure() 
ax = fig.add_subplot() 
ax.plot(dates, means, lw=1, color=martaBlue)
ax.fill_between(dates, means-stds, means+std, solor=martaBlue, alpha=0.5)
ax.set_xlabel("Time")
ax.set_ylabel("Power (kW)")

plt.tight_layout() 
plt.show() 