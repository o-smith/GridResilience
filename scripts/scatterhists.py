import numpy as np 
import matplotlib 
import matplotlib.pyplot as plt 
from scipy.stats import gaussian_kde

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
Icolour = '#DB2420' #Central line red
Ocolour = '#00A0E2' #Victoria line blue
O2colour = '#868F98' #Jubilee line grey
D2colour = '#F386A0' #Hammersmith line pink 
D4colour = '#97015E' #Metropolitan line magenta
D6colour = '#B05F0F' #Bakerloo line brown
D5colour = '#00843D' #District line green

# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
spacing = 0.005


rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]

fnames = [] 
fnames.append("data/cascade_s_vs_t_q75.txt")
fnames.append("data/cascade_s_vs_t_q10.txt")
fnames.append("data/cascade_s_vs_t_q121.txt")
fnames.append("data/cascade_s_vs_t_q20.txt") 
alphavals = [0.75, 1, 1.25, 2]

for i, fname in enumerate(fnames):

    # start with a rectangular Figure
    plt.figure()

    ax_scatter = plt.axes(rect_scatter)
    ax_scatter.tick_params(direction='in', top=True, right=True)
    ax_histx = plt.axes(rect_histx)
    ax_histx.tick_params(direction='in', labelbottom=False)
    ax_histy = plt.axes(rect_histy)
    ax_histy.tick_params(direction='in', labelleft=False)

    dat = np.genfromtxt(fname)
    y, x = np.hsplit(dat,2)
    x = x[:,0]
    y = y[:,0]
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    # the scatter plot:
    ax_scatter.scatter(x, y, c=z, s=30) #, edgecolor='', cmap="viridis")

    # now determine nice limits by hand:
    binwidth = 0.25
    lim = np.ceil(np.abs([x, y]).max() / binwidth) * binwidth
    ax_scatter.set_xlim((-1, 101))
    ax_scatter.set_ylim((-0.05, 1.05))

    ax_histx.hist(x, bins=70, density=True, width=max(x)/50.0, color=martaRed)
    ax_histy.hist(y, bins=70, orientation='horizontal', density=True, height=max(y)/50,color=martaBlue)

    ax_histx.set_xlim(ax_scatter.get_xlim())
    ax_histy.set_ylim(ax_scatter.get_ylim())

    plt.ylabel("$\\mathcal{S}$", rotation=0, labelpad=340)
    plt.xlabel("$T$")
    title = "$\\alpha=$" + str(alphavals[i])
    plt.title(title)
    plt.tight_layout() 
    plt.show()



