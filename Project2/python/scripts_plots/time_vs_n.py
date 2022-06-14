import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit
import sys

input_file = sys.argv[1]
output_image = sys.argv[2]
input_data = np.loadtxt(input_file)

x = input_data[1:,1]
y = input_data[1:,0]

x = np.log(x)
y = np.log(y)
slope, intercept, r, p, std_err = stats.linregress(x,y)
print (slope)

plt.plot(x,y)
plt.xlabel("log(K)",fontsize="20")
plt.ylabel("log(time)",fontsize="20")
plt.savefig(output_image)
plt.show()
