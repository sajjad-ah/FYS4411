import numpy as np
import matplotlib.pyplot as plt
import sys

nameFile = sys.argv[1]

dataFile = np.loadtxt(str(nameFile))
x = dataFile[:,0]
y = dataFile[:,1]
const = dataFile[:,2]

c = const[0]

plt.plot(x,y,"-o")
plt.plot(x,const,":")
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$<E_{L}>$")
plt.ylim(bottom = c-1, top = c+10)
plt.legend(["Average energy","Exact energy"])
plt.show()
