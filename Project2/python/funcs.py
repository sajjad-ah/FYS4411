#export OMP_NUM_THREADS=4

import numpy as np
import time
import sparse
import numba

@numba.njit(cache=True)
def remove_small_elements():
    return 0

def create_e_abij(occ, vir, e):
    e_abij = np.zeros((vir,vir,occ,occ))
    for a in range(occ,occ+vir):
        for b in range(occ,occ+vir):
            for i in range(occ):
                for j in range(occ):
                    e_abij[a-occ][b-occ][i][j] = float(1)/(e[i] + e[j] - e[a] - e[b])
    return e_abij

def transform_one_body_elements_sparse(h, c):
    _h = np.einsum("ip, jq, ij -> pq", c.conj(), c, h)

    _h[np.abs(_h) < 1e-8] = 0

    return sparse.COO.from_numpy(_h)

def transform_two_body_elements_sparse(u, c):
    _u = np.einsum("ls, ijkl -> ijks", c, u)
    _u = np.einsum("kr, ijks -> ijrs", c, _u)
    _u = np.einsum("jq, ijrs -> iqrs", c.conj(), _u)
    _u = np.einsum("ip, iqrs -> pqrs", c.conj(), _u)

    _u[np.abs(_u) < 1e-8] = 0

    return sparse.COO.from_numpy(_u)

def transform_two_body_elements(u, c):
    _u = np.einsum("ls, ijkl -> ijks", c, u)
    _u = np.einsum("kr, ijks -> ijrs", c, _u)
    _u = np.einsum("jq, ijrs -> iqrs", c.conj(), _u)
    _u = np.einsum("ip, iqrs -> pqrs", c.conj(), _u)
    return _u

def transform_2b(u,C):
    v = np.tensordot(C, u, axes=[(0),(3)])
    v = np.tensordot(C, v, axes=[(0),(3)])
    v = np.tensordot(C, v, axes=[(0),(3)])
    v = np.tensordot(C, v, axes=[(0),(3)])
    return v

def transform_2b_sparse(u,C):
    print (u.nbytes)
    u = sparse.tensordot(C, u, axes=[(0),(3)])
    print (u.nbytes)
    print (np.amin(abs(u.data)))
    u = sparse.tensordot(C, u, axes=[(0),(3)])
    u = sparse.tensordot(C, u, axes=[(0),(3)])
    u = sparse.tensordot(C, u, axes=[(0),(3)])
    print("2")
    return u

def transform_one_body_elements(h, c):
    _h = np.einsum("ip, jq, ij -> pq", c.conj(), c, h)
    return _h



def HF_transformed(int_1b,int_2b,occ):
    energy_1b = 0
    energy_2b = 0

    for i in range(occ):
        energy_1b += int_1b[i][i]

    for i in range(occ):
        for j in range(occ):
            energy_2b += int_2b[i][j][i][j]

    HF_energy = energy_1b + 0.5*energy_2b
    return HF_energy


##### correlation energy calculations #####
#@numba.njit(cache=True)

def corr_energy(occ,vir,t,int_2b):
    slice_occ = slice(0,occ)
    slice_vir = slice(occ,vir+occ)
    energy = 0.25*np.tensordot(int_2b[slice_occ,slice_occ,slice_vir,slice_vir],t,axes=[(0,1,2,3),(2,3,0,1)])
    return energy

def corr_energy_sparse(occ,vir,t,int_2b):
    slice_occ = slice(0,occ)
    slice_vir = slice(occ,vir+occ)
    energy = 0.25*sparse.tensordot(int_2b[slice_occ,slice_occ,slice_vir,slice_vir],t,axes=[(0,1,2,3),(2,3,0,1)])
    energy = energy.todense()
    return energy


def MP2(int_2b, e_abij, occ, vir):
    slice_occ = slice(0,occ)
    slice_vir = slice(occ,vir+occ)
    t = np.multiply(int_2b[slice_vir,slice_vir,slice_occ,slice_occ],e_abij)
    MP2_energy = corr_energy(occ,vir,t,int_2b)
    return t, MP2_energy

def MP2_sparse(int_2b, e_abij, occ, vir):
    slice_occ = slice(0,occ)
    slice_vir = slice(occ,vir+occ)
    t = np.multiply(int_2b[slice_vir,slice_vir,slice_occ,slice_occ],e_abij)
    print ("MP2_sparse",t.density)
    MP2_energy = corr_energy_sparse(occ,vir,t,int_2b)
    return t, MP2_energy

def CCD_solver_sparse(int_2b,e,occ,vir,t,E_ref,e_abij,opt=False):
    print ('int_2b',int_2b.density)
    print ('CCD_solver_1',t.density)
    t_ = 1*t
    print ('CCD_solver_2',t_.density)
    slice_occ = slice(0,occ)
    slice_vir = slice(occ,vir+occ)
    counter = 0
    eps = 10e-7
    dE = eps+1
    max_it = 100
    E_old = E_ref
    time_start = time.time()
    while (dE > eps) and ( counter < max_it ):
        print ("s1")
        #term1
        term = 1*int_2b[slice_vir,slice_vir,slice_occ,slice_occ]
        #term2
        #time_init = time.time()
        term += 0.5*sparse.tensordot(int_2b[slice_vir,slice_vir,slice_vir,slice_vir],t_,axes=[(2,3),(0,1)])
        print ("s2")
        #print (time.time()-time_init)
        #term3
        term += 0.5*sparse.tensordot(int_2b[slice_occ,slice_occ,slice_occ,slice_occ],t_,axes=[(0,1),(2,3)]).transpose((2,3,0,1))
        #term4
        temp = sparse.tensordot(int_2b[slice_occ,slice_vir,slice_vir,slice_occ],t_,axes=[(0,2),(3,1)]).transpose((2,0,3,1))
        term += temp - temp.transpose((0,1,3,2)) - temp.transpose((1,0,2,3)) + temp.transpose((1,0,3,2))
        #term
        inter = sparse.tensordot(int_2b[slice_occ,slice_occ,slice_vir,slice_vir],t_,axes=[(2,3),(0,1)])
        term += 0.25*sparse.tensordot(t_,inter,axes=[(2,3),(0,1)])
        #term6
        inter = sparse.tensordot(int_2b[slice_occ,slice_occ,slice_vir,slice_vir],t_,axes=[(0,2),(3,1)])
        temp = sparse.tensordot(inter,t_,axes=[(0,1),(3,1)]).transpose((0,2,1,3))
        term += temp - temp.transpose((0,1,3,2))
        #term7
        inter = sparse.tensordot(int_2b[slice_occ,slice_occ,slice_vir,slice_vir],t_,axes=[(0,2,3),(3,1,0)])
        temp = -0.5*sparse.tensordot(inter,t_,axes=[(0),(2)]).transpose((1,2,0,3))
        term += temp - temp.transpose((0,1,3,2))
        #term8
        inter = sparse.tensordot(int_2b[slice_occ,slice_occ,slice_vir,slice_vir],t_,axes=[(0,1,2),(3,2,1)])
        temp = -0.5*sparse.tensordot(inter,t_,axes=[(0),(0)])
        term += temp - temp.transpose((1,0,2,3))
        ##########
        #####print (type(term))
        print ('CCD_solver_3',t_.density)
        t_ = sparse.elemwise(np.multiply, e_abij,term)
        print ('CCD_solver_4',t_.density)
        CCD_energy = corr_energy_sparse(occ,vir,t_,int_2b)
        dE = abs(CCD_energy - E_old)
        E_old = CCD_energy
        counter += 1
    time_elapsed = (time.time() - time_start)
    print ("Time elapsed (CCD solver):  ",time_elapsed)
    print ("Number of iterations:       ",counter)
    return t_,CCD_energy

def CCD_solver_ten(int_2b,e,occ,vir,t,E_ref,e_abij,opt=False):
    t_old = np.copy(t)
    t_new = np.copy(t)
    slice_occ = slice(0,occ)
    slice_vir = slice(occ,vir+occ)
    counter = 0
    time_start = time.time()
    eps = 10e-7
    dE = eps+1
    max_it = 100
    E_old = E_ref
    while (dE > eps) and ( counter < max_it ):
        #term1
        term = np.copy(int_2b[slice_vir,slice_vir,slice_occ,slice_occ])
        #term2
        term += 0.5*np.tensordot(int_2b[slice_vir,slice_vir,slice_vir,slice_vir],t_old,axes=[(2,3),(0,1)])
        #term3
        term += 0.5*np.tensordot(int_2b[slice_occ,slice_occ,slice_occ,slice_occ],t_old,axes=[(0,1),(2,3)]).transpose(2,3,0,1)
        #term4
        temp = np.tensordot(int_2b[slice_occ,slice_vir,slice_vir,slice_occ],t_old,axes=[(0,2),(3,1)]).transpose(2,0,3,1)
        term += temp - temp.transpose(0,1,3,2) - temp.transpose(1,0,2,3) + temp.transpose(1,0,3,2)
        #term5
        inter = np.tensordot(int_2b[slice_occ,slice_occ,slice_vir,slice_vir],t_old,axes=[(2,3),(0,1)])
        term += 0.25*np.tensordot(t_old,inter,axes=[(2,3),(0,1)])
        #term6
        inter = np.tensordot(int_2b[slice_occ,slice_occ,slice_vir,slice_vir],t_old,axes=[(0,2),(3,1)])
        temp = np.tensordot(inter,t_old,axes=[(0,1),(3,1)]).transpose(0,2,1,3)
        term += temp - temp.transpose(0,1,3,2)
        #term7
        inter = np.tensordot(int_2b[slice_occ,slice_occ,slice_vir,slice_vir],t_old,axes=[(0,2,3),(3,1,0)])
        temp = -0.5*np.tensordot(inter,t_old,axes=[(0),(2)]).transpose(1,2,0,3)
        term += temp - temp.transpose(0,1,3,2)
        #term8
        inter = np.tensordot(int_2b[slice_occ,slice_occ,slice_vir,slice_vir],t_old,axes=[(0,1,2),(3,2,1)])
        temp = -0.5*np.tensordot(inter,t_old,axes=[(0),(0)])
        term += temp - temp.transpose(1,0,2,3)
        ##########
        t_new = np.multiply(e_abij,term)
        CCD_energy = corr_energy(occ,vir,t_new,int_2b)
        dE = abs(CCD_energy - E_old)
        E_old = CCD_energy
        t_old = np.copy(t_new)
        counter += 1
    time_elapsed = (time.time() - time_start)
    print ("Time elapsed (CCD solver):  ",time_elapsed)
    print ("Number of iterations:       ",counter)
    return t_new,CCD_energy

def CCD_solver_ein(int_2b,e,occ,vir,t,E_ref,e_abij,opt=False):
    theta = 0.5
    t_old = np.copy(t)
    t_new = np.copy(t)
    slice_occ = slice(0,occ)
    slice_vir = slice(occ,vir+occ)
    counter = 0
    eps = 10e-7
    dE = eps+1
    max_it = 1
    E_old = E_ref
    time_start = time.time()
    while ( counter < max_it ) and (dE > eps):
        #term1
        term = np.copy(int_2b[slice_vir,slice_vir,slice_occ,slice_occ])
        #term2
        term += 0.5*np.einsum("abcd, cdij -> abij",int_2b[slice_vir,slice_vir,slice_vir,slice_vir],t_old)
        #term3
        term += 0.5*np.einsum("klij, abkl -> abij",int_2b[slice_occ,slice_occ,slice_occ,slice_occ],t_old)
        #term4
        temp = np.einsum("kbcj, acik -> abij",int_2b[slice_occ,slice_vir,slice_vir,slice_occ],t_old)
        term += temp - temp.transpose(0,1,3,2) - temp.transpose(1,0,2,3) + temp.transpose(1,0,3,2)
        #term5
        term += 0.25*np.einsum("klcd, cdij, abkl -> abij",int_2b[slice_occ,slice_occ,slice_vir,slice_vir],t_old,t_old,optimize=opt)
        #term6
        temp = np.einsum("klcd, acik, bdjl -> abij",int_2b[slice_occ,slice_occ,slice_vir,slice_vir],t_old,t_old,optimize=opt)
        term += temp - temp.transpose(0,1,3,2)
        #term7
        temp = -0.5*np.einsum("klcd, dcik, ablj -> abij",int_2b[slice_occ,slice_occ,slice_vir,slice_vir],t_old,t_old,optimize=opt)
        term += temp - temp.transpose(0,1,3,2)
        #term8
        temp = -0.5*np.einsum("klcd, aclk, dbij -> abij",int_2b[slice_occ,slice_occ,slice_vir,slice_vir],t_old,t_old,optimize=opt)
        term += temp - temp.transpose(1,0,2,3)
        t_new = np.multiply(e_abij,term)

        CCD_energy = corr_energy(occ,vir,t_new,int_2b)
        dE = abs(CCD_energy - E_old)
        E_old = CCD_energy
        t_old = np.copy(t_new)
        counter += 1
    time_elapsed = (time.time() - time_start)
    print ("Time elapsed (CCD solver):  ",time_elapsed)
    print ("Number of iterations:       ",counter)
    return t_new,CCD_energy




"""
Time elapsed (CCD solver):   0.07738600000000001
Number of iterations:        10
Total CCD energy (EINSUM):        (3.039047993881479+0j)
CCD correlation energy (EINSUM):  (-0.1236433559841587+0j)
"""


##### Old, non-used functions #####

def CCD_solver(int_2b,e,occ,vir,t,E_ref):
    t_old = np.copy(t)
    t_new = np.copy(t)
    counter = 0
    time_start = time.time()
    eps = 10e-7
    dE = eps+1
    max_it = 10
    E_old = E_ref
    while (dE > eps) and ( counter < max_it ):
        for a in range(occ,occ+vir):
            for b in range(occ,occ+vir):
                for i in range(occ):
                    for j in range(occ):
                        ### Term 1 ###
                        term_1 = int_2b[a][b][i][j]
                        ### Term 2 ###
                        term_2 = 0
                        for c in range(occ,occ+vir):
                            for d in range(occ,occ+vir):
                                term_2 += int_2b[a][b][c][d]*t_old[c-occ][d-occ][i][j]
                        term_2 = 0.5*term_2
                        ### Term 3 ###
                        term_3 = 0
                        for k in range(occ):
                            for l in range(occ):
                                term_3 += int_2b[k][l][i][j]*t_old[a-occ][b-occ][k][l]
                        term_3 = 0.5*term_3
                        ### Term 4 ###
                        term_4 = 0
                        for k in range(occ):
                            for c in range(occ,occ+vir):
                                term_4 += int_2b[k][b][c][j]*t_old[a-occ][c-occ][i][k] \
                                        - int_2b[k][b][c][i]*t_old[a-occ][c-occ][j][k] \
                                        - int_2b[k][a][c][j]*t_old[b-occ][c-occ][i][k] \
                                        + int_2b[k][a][c][i]*t_old[b-occ][c-occ][j][k]
                        ### Term 5, 6, 7 and 8 ###
                        term_5 = 0
                        term_6 = 0
                        term_7 = 0
                        term_8 = 0
                        for k in range(occ):
                            for l in range(occ):
                                for c in range(occ,occ+vir):
                                    for d in range(occ,occ+vir):
                                        term_5 += int_2b[k][l][c][d]*t_old[c-occ][d-occ][i][j]*t_old[a-occ][b-occ][k][l]
                                        term_6 += int_2b[k][l][c][d] \
                                        *(t_old[a-occ][c-occ][i][k]*t_old[b-occ][c-occ][j][l] - t_old[a-occ][c-occ][j][k]*t_old[b-occ][c-occ][i][l])
                                        term_7 += int_2b[k][l][c][d] \
                                        *(t_old[d-occ][c-occ][i][k]*t_old[a-occ][b-occ][l][j] - t_old[d-occ][c-occ][j][k]*t_old[a-occ][b-occ][l][i])
                                        term_8 += int_2b[k][l][c][d] \
                                        *(t_old[a-occ][c-occ][l][k]*t_old[d-occ][b-occ][i][j] - t_old[b-occ][c-occ][l][k]*t_old[d-occ][a-occ][i][j])
                        term_5 = 0.25*term_5
                        term_7 = -0.5*term_7
                        term_8 = -0.5*term_8
                        tot_term = term_1 + term_2 + term_3 + term_4 \
                        + term_5 + term_6 + term_7 + term_8
                        t_new[a-occ][b-occ][i][j] = tot_term / float(e[i]+e[j]-e[a]-e[b])
        CCD_energy = corr_energy(occ,vir,t_new,int_2b)
        dE = abs(CCD_energy - E_old)
        E_old = CCD_energy
        t_old = np.copy(t_new)
        counter += 1
    time_elapsed = (time.time() - time_start)
    print ("Time elapsed (CCD solver):  ",time_elapsed)
    print ("Number of iterations:       ",counter)
    return t_new, CCD_energy


def primitive_2b_transform(u, C):
    num_basis = len(C[0,:])
    u_ = np.zeros((num_basis,num_basis,num_basis,num_basis),dtype=complex)
    i = 0
    for p in range(num_basis):
        i+=1
        percent = float(i-1)/num_basis*100
        print ("%f %%" %(percent))
        for q in range(num_basis):
            for r in range(num_basis):
                for s in range(num_basis):
                    sum = 0
                    for alpha in range(num_basis):
                        for beta in range(num_basis):
                            for gamma in range(num_basis):
                                for delta in range(num_basis):
                                    sum += C[alpha][p]*C[beta][q]*C[gamma][r]*C[delta][s]*u[alpha][beta][gamma][delta]
                    u_[p][q][r][s] = sum
    return u_

def MP2_old(int_2b, e, occ, vir):
    t = np.zeros((vir,vir,occ,occ),dtype=complex)
    for a in range(occ,occ+vir):
        for b in range(occ,occ+vir):
            for i in range(occ):
                for j in range(occ):
                    t[a-occ][b-occ][i][j] = int_2b[a][b][i][j] / (e[i]+e[j]-e[a]-e[b])
    MP2_energy = corr_energy(occ,vir,t,int_2b)
    return t, MP2_energy

def MP2_sparse_old(int_2b, e, occ, vir):
    t = np.zeros((vir,vir,occ,occ),dtype=complex)
    for a in range(occ,occ+vir):
        for b in range(occ,occ+vir):
            for i in range(occ):
                for j in range(occ):
                    t[a-occ][b-occ][i][j] = int_2b[a][b][i][j] / (e[i]+e[j]-e[a]-e[b])
    MP2_energy = corr_energy_sparse(occ,vir,t,int_2b)
    t[np.abs(t) < 1e-8] = 0
    return sparse.COO.from_numpy(t), MP2_energy

def corr_energy_old(occ,vir,t,int_2b):
    #print (type(occ))
    energy = 0
    for a in range(occ,occ+vir):
        for b in range(occ,occ+vir):
            for i in range(occ):
                for j in range(occ):
                    energy += int_2b[i][j][a][b]*t[a-occ][b-occ][i][j]
    energy = 0.25*energy
    return energy
