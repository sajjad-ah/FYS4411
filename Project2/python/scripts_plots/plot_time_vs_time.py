import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import stats
from scipy.optimize import curve_fit

input_data = np.loadtxt(sys.argv[1])
ein = input_data[:,0]
ten = input_data[:,1]
num_basis = []
ein_data = []
ten_data = []

c = 0
for i in range(int(len(input_data)/3.0)):
    val_ein = 0
    val_ten = 0
    for j in range(3):
        val_ein += ein[c+j]
        val_ten += ten[c+j]
    val_ein = val_ein/3.0
    val_ten = val_ten/3.0
    val_ein = val_ein/float(input_data[c,3])
    val_ten = val_ten/float(input_data[c,3])
    num_basis.append(input_data[c,2])
    ein_data.append(val_ein)
    ten_data.append(val_ten)
    c += 3

input_data2 = np.loadtxt(sys.argv[2])
ein_false = input_data2[:,0]
num_basis_false = []
ein_data_false = []

c = 0
for i in range(int(len(input_data2)/3.0)):
    val_ein = 0
    for j in range(3):
        val_ein += ein_false[c+j]
    val_ein = val_ein/3.0
    val_ten = val_ten/3.0
    val_ein = val_ein/float(input_data2[c,2])
    num_basis_false.append(input_data2[c,1])
    ein_data_false.append(val_ein)
    c += 3

ein_data = np.array(ein_data)
ten_data = np.array(ten_data)
ein_data_false = np.array(ein_data_false)


ein_data = np.log(ein_data)
ten_data = np.log(ten_data)
ein_data_false = np.log(ein_data_false)

num_basis = np.array(num_basis)
num_basis = np.log(num_basis)
num_basis_false = np.array(num_basis_false)
num_basis_false = np.log(num_basis_false)

slope1, intercept1, r1, p1, std_err1 = stats.linregress(num_basis[2:],ein_data[2:])
slope2, intercept2, r2, p2, std_err2 = stats.linregress(num_basis[2:],ten_data[2:])
slope3, intercept3, r3, p3, std_err3 = stats.linregress(num_basis_false[2:],ein_data_false[2:])

print ("Einslope: ", slope1)
print ("Tenslope: ", slope2)
print ("Einslope false: ", slope3)

plt.plot(num_basis,ein_data, '--o' ,label="Einsum")
plt.plot(num_basis,ten_data, '--o' ,label="Tensordot")
plt.plot(num_basis_false,ein_data_false, '--o' ,label="Einsum (optimize=False)")
plt.xlabel("log(Number of basis functions)")
plt.ylabel("log(time)")
plt.legend()
plt.savefig("ein_vs_ten.png")
plt.show()
