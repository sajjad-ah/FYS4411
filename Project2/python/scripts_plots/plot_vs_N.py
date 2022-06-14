import numpy as np
import matplotlib.pyplot as plt

opt_true = np.loadtxt("N2_optTrue.dat")
opt_false = np.loadtxt("N2_optFalse.dat")

y_true = opt_true[:,4]/opt_true[:,5]
y_false = opt_false[:,4]/opt_false[:,5]

y_true = np.log(y_true)
y_false = np.log(y_false)

x = opt_true[:,0]
x = np.log(x)

plt.plot(x,y_true,label='True')
plt.plot(x,y_false,label='False')
plt.legend()
plt.show()
