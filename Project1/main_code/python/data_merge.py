import numpy as np
import sys

a = np.loadtxt("0.dat")
b = np.loadtxt("1.dat")
c = np.loadtxt("2.dat")
d = np.loadtxt("3.dat")

e = np.concatenate((a,b),axis=0)
f = np.concatenate((c,d),axis=0)
g = np.concatenate((e,f),axis=0)

np.savetxt(sys.argv[1],g)
