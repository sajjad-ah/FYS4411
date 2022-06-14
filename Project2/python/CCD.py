import numpy as np
from funcs import *
import time
import sys

if len(sys.argv) == 1:
    folder = ""
else:
    folder = sys.argv[1]

theta = 0
theta_list = [0.1,0.2,0.3]

C = np.load(folder+"C_matrix.npy")
h = np.load(folder+"h_matrix.npy")
u = np.load(folder+"u_matrix.npy")
e = np.load(folder+"e.npy")
occ = np.load(folder+"primitive_occ.npy")
occ = int(occ)

num_basis = len(C[0,:])
vir =  num_basis - occ

e_abij = create_e_abij(occ,vir,e)


int_1b = transform_one_body_elements(h, C)

int_2b = transform_2b(u, C)


###computing HF energy###
HF_energy = HF_transformed(int_1b,int_2b,occ)

print()
print("##### HF #####")
print ("HF reference energy:    ",HF_energy)

###computing MP2 amplitudes and energy###
t, MP2_energy = MP2(int_2b,e_abij,occ,vir)

print (int_2b.nbytes)

print()
print("##### MP2 #####")
print ("Total MP2 energy:       ",MP2_energy+HF_energy)
print ("MP2 correlation energy: ",MP2_energy)

print()
print("CCD calculation has started...")

#print()
#print("##### CCD einsum #####")
#t_CCD_ein, CCD_energy_ein = CCD_solver_ein(int_2b,e,occ,vir,t,MP2_energy,e_abij,opt=True)
#print ("Total CCD energy (EINSUM):       ",CCD_energy_ein+HF_energy)
#print ("CCD correlation energy (EINSUM): ",CCD_energy_ein)

print()
print("##### CCD tensordot #####")
t_CCD_ten, CCD_energy_ten = CCD_solver_ten(int_2b,e,occ,vir,t,MP2_energy,e_abij,theta,opt=True)
print ("Theta:                              ",theta)
print ("Total CCD energy (TENSORDOT):       ",CCD_energy_ten+HF_energy)
print ("CCD correlation energy (TENSORDOT): ",CCD_energy_ten)

#comment out this if only theta=0 run
for i in range(len(theta_list)):
    print()
    print("##### CCD tensordot #####")
    t_CCD_ten, CCD_energy_ten = CCD_solver_ten(int_2b,e,occ,vir,t,MP2_energy,e_abij,theta_list[i],opt=True)
    print ("Theta:                              ",theta_list[i])
    print ("Total CCD energy (TENSORDOT):       ",CCD_energy_ten+HF_energy)
    print ("CCD correlation energy (TENSORDOT): ",CCD_energy_ten)

###computing CCD amplitudes and energy###
#print()
#print("##### CCD #####)")
#t_CCD, CCD_energy = CCD_solver(int_2b,e,occ,vir,t,MP2_energy)
#print ("Total CCD energy:       ",CCD_energy+HF_energy)
#print ("CCD correlation energy: ",CCD_energy)
