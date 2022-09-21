from eigentools import Eigenproblem
from dedalus import public as de
import numpy as np
import scipy
import logging
logger = logging.getLogger(__name__)


def LOP_setup_mats(EP, parameters, norm, pencil=0, dense=True, k=100):
    if dense:
        EP.solve(sparse=False, pencil=pencil, parameters=parameters)
    else:
        EP.solve(sparse=True, N=k, pencil=pencil, parameters=parameters)
    if len(EP.evalues_good) == 0:
        raise ValueError("No good eigenvalues for these parameters")
    pre_right = EP.solver.pencils[pencil].pre_right  # right-preconditioner for the matrix pencil, call it PR
    pre_right_LU = scipy.sparse.linalg.splu(pre_right.tocsc())  # LU factorization of the right-preconditioner

    # eigenvalues = EP.evalues
    # if dense:
        # eigenvectors = ((EP.solver.eigenvectors.T)[EP.evalues_index]).T  # pretty sure there's a better way to do this
    # else:
        # eigenvectors = EP.solver.eigenvectors
    eigenvalues = EP.evalues_good
    eigenvectors = ((EP.solver.eigenvectors.T)[EP.evalues_index]).T
    V = pre_right_LU.solve(eigenvectors)  # V satisfies PR @ V = solver.eigenvectors. I think this means V is the eigenvectors WITHOUT pre-conditioning
    Q, R = np.linalg.qr(V)  # Q is an orthonormal (in terms of the 2-norm!) basis in the subspace spanned by the eigenvectors in V, R is the change of basis tsfm
    M_E = EP.compute_mass_matrix(pre_right @ Q, norm)  # This is M represented in the basis of the columns of Q
    Fdagger = np.linalg.cholesky(M_E)
    F = Fdagger.conj().T  # this gives ||F @ u||_2 = ||u||_E for u expressed in the (k-dimensional!) basis given by Q

    return eigenvalues, eigenvectors, F, Q, R, pre_right


def LOP(EP, times, parameters, norm, pencil=0, dense=True, k=100):

    Gs = np.zeros_like(times)
    eigenvalues, eigenvectors, F, Q, R, pre_right = LOP_setup_mats(EP, parameters, norm, pencil=pencil, dense=dense, k=k)
    for tind, t in enumerate(times):
        K = R @ np.diag(np.exp(eigenvalues * t)) @ np.linalg.inv(R)
        s_vals = scipy.linalg.svdvals(F @ K @ np.linalg.inv(F))
        Gs[tind] = s_vals[0]
    Gs = Gs**2.0

    return Gs, eigenvalues, eigenvectors


def F_from_norm(EP, pre_right, Q, norm, ma):  # I want this to be its own function so I can calculate F_KE, F_ME, etc.
    if EP.evalues is None:
        raise NotImplementedError("F_from_norm assumes EP is already solved.")
    M_E = EP.compute_mass_matrix(pre_right @ Q, lambda x1, x2: norm(x1, x2, ma))  # when norm is KE, ma is totally unnecessary. How do I generalize this with **kwargs magic?
    Fdagger = np.linalg.cholesky(M_E)
    F = Fdagger.conj().T
    return F
