import warnings
import scipy.linalg as scl


class HartreeFock:
    def __init__(self, system, verbose=False, np=None):
        if np is None:
            import numpy as np

        self.np = np
        self.system = system
        self.verbose = verbose

    def scf(self, max_iters=100, tolerance=1e-9):
        np = self.np

        occ = self.system.o
        h, u = self.system.h, self.system.u
        ##########################################################
        # Treat all basis sets as non-orthogonal
        s, U = np.linalg.eigh(self.system.s)  # Eq. (3.166) Szabo-Ostlund
        X = np.dot(
            U * 1.0 / np.sqrt(s), np.conj(U).T
        )  # Eq. (3.167) Szabo-Ostlund

        # Calculate initial core guess: [Szabo:1996] pp. 145
        Hp = X @ h @ X  # Eqn. 3.177
        e, C2 = np.linalg.eigh(Hp)  # Solving Eqn. 1.178
        self.C = X @ C2  # Back transform, Eqn. 3.174
        ###########################################################

        # Reference energy as the initial energy
        e_hf_prev = np.einsum("ii->", h[occ, occ]) + 0.5 * np.einsum(
            "ijij->", u[occ, occ, occ, occ]
        )
        # SCF-procedure
        for i in range(max_iters + 1):

            D = np.dot(self.C[:, occ], np.conj(self.C[:, occ]).T)
            F = h + np.einsum("sr,psqr->pq", D, u)

            Fp = X.dot(F).dot(X)
            eps, C2 = np.linalg.eigh(Fp)
            self.C = X.dot(C2)

            self.e_hf = self._hf_energy(D)

            if abs(self.e_hf - e_hf_prev) < tolerance:
                break

            if i < max_iters:
                e_hf_prev = self.e_hf

        converged = abs(self.e_hf - e_hf_prev) < tolerance
        if converged and self.verbose:
            print("HF converged to given precision in {0} iterations".format(i))

        if not converged:
            warnings.warn(
                "HF did not converge to given precision {0}".format(tolerance)
            )

        if self.verbose:
            print("Ehf: {0}".format(self.e_hf))

        np.save("h_matrix",h)
        np.save("u_matrix",u)
        np.save("occ",occ)
        np.save("e",eps)
        return self.C

    def _hf_energy(self, D):
        np = self.np

        occ = self.system.o
        h, u = self.system.h, self.system.u

        D_new = np.dot(self.C[:, occ], np.conj(self.C[:, occ]).T)

        return np.einsum("pq,pq->", D, h) + 0.5 * np.einsum(
            "ag,bd,abgd->", D, D_new, u
        )


from quantum_systems import *
import numpy as np

#system setup
num_part = 2
num_basis = 6
tdho = TwoDimensionalHarmonicOscillator(num_part,num_basis,5,1,omega=1.0) #(num_of_particles,num_of_basis_functions,harmonic_oscillstor_size,grid_points)
tdho.setup_system()

#run HF SCF
HF1 = HartreeFock(tdho, verbose=True)


np.save("C_matrix",HF1.scf())
np.save("primitive_occ",num_part)
