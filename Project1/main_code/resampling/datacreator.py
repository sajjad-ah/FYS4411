import numpy as np
import random
import sys

num_of_data_points = int(sys.argv[1])

datapoints = np.zeros(num_of_data_points)

for i in range(len(datapoints)):
    datapoints[i] = 10*random.random()

np.savetxt("datapoints.dat",datapoints)
