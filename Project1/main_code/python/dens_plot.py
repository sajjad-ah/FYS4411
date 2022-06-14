import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return np.exp(-1.0*x**2)

a = np.loadtxt("dens_noint.dat")
b = np.loadtxt("dens_int1.dat")
c = np.loadtxt("dens_int20_second.dat")

w = np.zeros((len(c[:,0]),2))

for i in range(len(c[:,0])):
    w[i,0] = c[i,0]

for i in range(len(w[:,0])):
    w[i,1] = f(w[i,0])


a_sum = 0
b_sum = 0
c_sum = 0
w_sum = 0
for i in range(len(a[:,0])):
    a_sum += a[i,1]
    b_sum += b[i,1]
    c_sum += c[i,1]
    w_sum += w[i,1]

a[:,1] = a[:,1]/a_sum
b[:,1] = b[:,1]/b_sum
c[:,1] = c[:,1]/c_sum
w[:,1] = w[:,1]/w_sum

plt.plot(a[:,0],a[:,1],label="MC, non-interacting")
plt.plot(b[:,0],b[:,1],'g',label="MC, interacting")
plt.plot(c[:,0],c[:,1],label="int20")
plt.plot(w[:,0],w[:,1],'r--',label="Analytical, non-interacting")
plt.xlabel(r"Distance ($r/a_{ho}$)")
plt.ylabel(r"One-body density")
plt.legend()
plt.savefig("1body_dens.png")
plt.show()
