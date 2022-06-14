import numpy as np
import sys

a = np.loadtxt("0.pos")
b = np.loadtxt("1.pos")
c = np.loadtxt("2.pos")
d = np.loadtxt("3.pos")

e = np.concatenate((a,b),axis=0)
f = np.concatenate((c,d),axis=0)
g = np.concatenate((e,f),axis=0)

np.savetxt(sys.argv[1],g)
