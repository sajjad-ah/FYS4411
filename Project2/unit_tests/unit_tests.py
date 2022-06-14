import numpy as np
import time
import sys
sys.path.insert(0, '../python')
from funcs import *

C = np.load("UT_C_matrix.npy")
h = np.load("UT_h_matrix.npy")
u = np.load("UT_u_matrix.npy")
e = np.load("UT_e.npy")
occ = np.load("UT_primitive_occ.npy")
occ = int(occ)

num_basis = len(C[0,:])
vir =  num_basis - occ

e_abij = create_e_abij(occ,vir,e)

int_1b = transform_one_body_elements(h, C)

int_2b = transform_2b(u, C)

###
energy_from_HF = 3.2533141373155003+0j

###computing HF energy###
HF_energy = HF_transformed(int_1b,int_2b,occ)

if energy_from_HF == HF_energy:
    print ("HF energy:      PASSED")
else:
    print ("HF energy:      FAILED")
    print ("HF energy:      The energy in HO basis and HF basis are different")
    print ("HF energy:      The energy in HO basis: ",energy_from_HF)
    print ("HF energy:      The energy in HF basis: ",HF_energy)
