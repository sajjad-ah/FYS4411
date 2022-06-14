import numpy as np
from funcs_sparse import *
import time
import sys

eps = 1e-8
print ("1")
if len(sys.argv) == 1:
    folder = ""
else:
    folder = sys.argv[1]

h = np.load(folder+"h_matrix.npy")
C_data = np.load(folder+"C_data.npy")
C_coords = np.load(folder+"C_coords.npy")
u_data = np.load(folder+"u_data.npy")
u_coords = np.load(folder+"u_coords.npy")
e = np.load(folder+"e.npy")
occ = np.load(folder+"primitive_occ.npy")
num_basis = np.load(folder+"num_basis.npy")
occ = int(occ)
num_basis = int(num_basis)
vir =  num_basis - occ
C = sparse.COO(C_coords, C_data, shape=(occ+vir, occ+vir))
u = sparse.COO(u_coords, u_data, shape=(occ+vir,occ+vir,occ+vir, occ+vir))
print ("2")
#C.data = np.delete(C.data,[6])
#C.coords = np.delete(C.coords,6,1)

e_abij = create_e_abij(occ,vir,e)
print ("3")
int_1b_sparse = transform_one_body_elements_sparse(h, C)
print ("4")
#C[np.abs(C) < eps] = 0
#C = sparse.COO.from_numpy(C)
#int_2b_sparse = transform_2b_sparse(u, C)
int_2b_sparse = u
print (int_2b_sparse.nbytes)
int_2b_sparse.coords, int_2b_sparse.data = r(int_2b_sparse.coords,int_2b_sparse.data)
print (int_2b_sparse.nbytes)
print (int_2b_sparse.density)
#int_2b_sparse = u

###computing HF energy###
HF_energy = HF_transformed(int_1b_sparse,int_2b_sparse,occ)
HF_energy =
print()
print("##### HF #####")
print ("HF reference energy:    ",HF_energy)

###computing MP2 amplitudes and energy###
t_sparse, MP2_energy_sparse = MP2_sparse(int_2b_sparse,e_abij,occ,vir)


#t_sparse, MP2_energy_sparse = MP2_sparse_old(int_2b_sparse,e,occ,vir)


print()
print("##### MP2 #####")
print ("Total MP2 energy:       ",MP2_energy_sparse+HF_energy)
print ("MP2 correlation energy: ",MP2_energy_sparse)
#print ("MP2 correlation energy sparse: ",MP2_energy_sparse)

print()
print("CCD calculation has started...")

print()
print("##### CCD sparse #####)")
t_CCD_ten_sparse, CCD_energy_ten_sparse = CCD_solver_sparse(int_2b_sparse,e,occ,vir,t_sparse,MP2_energy_sparse,e_abij,opt=False)
print ("Total CCD energy (SPARSE):       ",CCD_energy_ten_sparse+HF_energy)
print ("CCD correlation energy (SPARSE): ",CCD_energy_ten_sparse)
