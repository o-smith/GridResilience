import matplotlib.pyplot as plt 
import numpy as np 
import matplotlib
from matplotlib import colors

#Set up plot style
font = {'size'   : 11}
matplotlib.rc('font', **font)
matplotlib.rc('font', serif='Computer Modern Roman')

z = 100
n = 200
alphamin = 0.01
alphamax = 3.0
alphares = 600
d = np.genfromtxt("data/rho_dostancek6_200_normed.txt")
d2 = np.genfromtxt("data/rho_dostancek10_200_normed.txt")

x = np.linspace(0.665, 3.0, 1000)
y = np.linspace(0.665, 3.0, 500)
w = np.linspace(0.665, 1.0, 500)
# print(x)
curve = (n-2)*x/((2+x))
curve2 = (n-2)*((1-w)/(2-w)) 
curve3 = (n-2)*(y-1)/y

fig = plt.figure(figsize=(7,3)) 
ax = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

im = ax.imshow(d, interpolation="none", cmap="jet", aspect="auto", extent=[alphamin,alphamax,0,z-2])
# ax.plot(x, curve, "w", lw=2, ls="--")
# ax.plot(w, curve2, "w", lw=2, ls="--")
# ax.plot(y, curve3, "w", lw=2, ls="--")
ax.axvline(x=6/5, color="w", ls="--", lw=2)
ax.set_xlim([alphamin,alphamax])
ax.set_ylim([0,z-2])
ax.set_xlabel("$\\alpha/\\alpha_{\\ast}$", fontsize=14)
ax.set_ylabel("$r$", fontsize=14, rotation=0, labelpad=20)
ax.text(-0.5,110, "(a)", fontsize=11)
plt.colorbar(im, ax=ax)

z = 100
n = 200

im = ax2.imshow(d2, interpolation="none", cmap="jet", aspect="auto", extent=[alphamin,alphamax,0,z-2])
# ax.plot(x, curve, "w", lw=2, ls="--")
# ax.plot(w, curve2, "w", lw=2, ls="--")
# ax.plot(y, curve3, "w", lw=2, ls="--")
ax2.axvline(x=10/9, color="w", ls="--", lw=2)
ax2.set_xlim([alphamin,alphamax])
ax2.set_ylim([0,z-2])
ax2.set_xlabel("$\\alpha/\\alpha_{\\ast}$", fontsize=14)
ax2.set_ylabel("$r$", fontsize=14, rotation=0, labelpad=20)
ax2.text(-0.5,110, "(b)", fontsize=11)
# ax.text(1.0,75, "(i)", color="w")
# ax.text(0.55,28, "(ii)", color="w")
# ax.text(1.22,22, "(iii)", color="w")
# ax.text(2.05,55, "(iv)", color="w")
fig.colorbar(im)


plt.tight_layout()  
plt.show() 





