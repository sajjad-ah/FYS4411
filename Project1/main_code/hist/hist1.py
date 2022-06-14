import numpy as np
import matplotlib.pyplot as plt
import sys

#
#N_points = 10000
#np.random.seed(6789)

bins = 100

#data1 = np.random.randn(N_points)

data1 = np.loadtxt(sys.argv[1])

#result = plt.hist(data1, bins=20, color='b', edgecolor='k', alpha=0.55)
result = plt.hist(data1, bins)
plt.axvline(data1.mean(), color='r', linestyle='dashed', linewidth=1)

_, max_ = plt.ylim()
plt.text(data1.mean() + data1.mean()/10,
         max_ - max_/10,
         'Mean r: {:.6f}'.format(data1.mean()))

plt.xlabel('$r/a_{ho}$')
plt.ylabel('Relative frequency')
plt.title(r'D =,  N = , MC cycles = , $\alpha$ = ')
plt.grid(True)

plt.show()
