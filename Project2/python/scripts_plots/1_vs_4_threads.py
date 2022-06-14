import numpy as np
import matplotlib.pyplot as plt

in_data = np.loadtxt("1_vs_4_threads.data")

plt.plot(in_data[:,0],in_data[:,1],label="Threads: 1")
plt.plot(in_data[:,0],in_data[:,2],label="Threads: 4")
plt.xlabel("Number of basis functions",fontsize="16")
plt.ylabel("Time per iteration (s)",fontsize="16")
plt.legend(fontsize="16")
plt.savefig("1_vs_4_threads.png")
plt.show()
