import matplotlib 
import matplotlib.pyplot as plt 
import numpy as np 

#Set up plot style
font = {'size'   : 22}
matplotlib.rc('font', **font)
matplotlib.rc('font', serif='Computer Modern Roman')

#Define colours 
martaRed = "#%02x%02x%02x" %(194, 76, 81)
martaGreen = "#%02x%02x%02x" %(84, 166, 102)
martaBlue = "#%02x%02x%02x" %(76, 112, 176)
martaPurple = "#%02x%02x%02x" %(127, 112, 176) 
martaGold = "#%02x%02x%02x" %(204, 184, 115) 


#Colours
Icolour = '#DB2420' #Central line red
Ocolour = '#00A0E2' #Victoria line blue
O2colour = '#868F98' #Jubilee line grey
D2colour = '#F386A0' #Hammersmith line pink
D4colour = '#97015E' #Metropolitan line magenta
D6colour = '#B05F0F' #Bakerloo line brown
D5colour = '#00843D' #District line green
O3colour = '#021EA9'


dat = np.genfromtxt("data/smq6_32_32_0_n1.txt")
n1, ac1, std1 = np.hsplit(dat, 3) 

dat = np.genfromtxt("data/smq6_32_32_0_n2.txt")
n2, ac2, std2 = np.hsplit(dat, 3) 

dat = np.genfromtxt("data/smq6_32_32_0_n3.txt")
n3, ac3, std3 = np.hsplit(dat, 3) 

# print(np.shape(n1[:,0]))
n = np.append(n1[:,0], n2[:,0])
ac = np.append(ac1[:,0], ac2[:,0])

n = np.append(n, n3[:,0])
ac = np.append(ac, ac3[:,0])

plt.plot(n, ac, ".-", color=martaBlue, lw=3.0, label="$q=0.6$") 


dat = np.genfromtxt("data/smq1_32_32_0_n1.txt")
n1, ac1, std1 = np.hsplit(dat, 3) 

dat = np.genfromtxt("data/smq1_32_32_0_n2.txt")
n2, ac2, std2 = np.hsplit(dat, 3) 


# print(np.shape(n1[:,0]))
n = np.append(n1[:,0], n2[:,0])
ac = np.append(ac1[:,0], ac2[:,0])

n = np.append(n, n3[:,0])
ac = np.append(ac, ac3[:,0])

plt.plot(n, ac, ".-", color=martaRed, lw=3.0, label="$q=0.1$") 
plt.xticks([10,150])
plt.yticks([0.04,0.16])
plt.xlabel("$n$", fontsize=34, labelpad=-20)
plt.ylabel("$\\overline{a_c}$", fontsize=34, labelpad=-20, rotation=0)
plt.legend() 
plt.tight_layout() 
plt.show() 




