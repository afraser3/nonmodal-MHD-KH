from eigentools import Eigenproblem
from dedalus import public as de
import numpy as np
import scipy
import logging
logger = logging.getLogger(__name__)

def LOP(EP, times, parameters, norm, pencil=0, dense=True):

    Gs = np.zeros_like(times)
    if dense:
        EP.solve(sparse=False, pencil=pencil, parameters=parameters)
    else:
        EP.solve(sparse=True, N=k, pencil=pencil, parameters=parameters)
    pre_right = EP.solver.pencils[pencil].pre_right  # right-preconditioner for the matrix pencil, call it PR
    pre_right_LU = scipy.sparse.linalg.splu(pre_right.tocsc())  # LU factorization of the right-preconditioner
    if dense:
        eigenvalues = EP.evalues_good
        eigenvectors = ((EP.solver.eigenvectors.T)[EP.evalues_index]).T  # pretty sure there's a better way to do this
    else:
        eigenvalues = EP.solver.eigenvalues
        eigenvectors = EP.solver.eigenvectors
    V = pre_right_LU.solve(eigenvectors)  # V satisfies PR @ V = solver.eigenvectors. I think this means V is the eigenvectors WITHOUT pre-conditioning
    Q, R = np.linalg.qr(V)  # Q is an orthonormal (in terms of the 2-norm!) basis in the subspace spanned by the eigenvectors in V, R is the change of basis tsfm
    M_E = EP.compute_mass_matrix(pre_right @ Q, lambda x1, x2: norm(x1, x2, ma))  # This is M represented in the basis of the columns of Q
    Fdagger = np.linalg.cholesky(M_E)
    F = Fdagger.conj().T  # this gives ||F @ u||_2 = ||u||_E for u expressed in the (k-dimensional!) basis given by Q
    for tind, t in enumerate(times):
        K = R @ np.diag(np.exp(eigenvalues * t)) @ np.linalg.inv(R)
        s_vals = scipy.linalg.svdvals(F @ K @ np.linalg.inv(F))
        Gs[tind] = s_vals[0]
    Gs = Gs**2.0

    return Gs, eigenvalues, eigenvectors
