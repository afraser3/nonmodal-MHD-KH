from eigentools import Eigenproblem
from dedalus import public as de
import numpy as np
import scipy
from matplotlib import pyplot as plt
import logging
logger = logging.getLogger(__name__)

Nz = 256
Reynolds = 50.0
Pm = 1.0
mReynolds = Pm * Reynolds
Lz = 10.0*np.pi
kx = 0.2
# MA = 10.0  # Alfven Mach number. For this value, the system is unstable.
MA = 1.0  # For this value, the system is stable.
MAs = np.array([0.5, 0.75])  # , 1.0, 1.25])  # , 1.5])
MA2 = MA**2.0
k = 250
dense = False
if dense:
    fname = 'LOP_MA_scan_Nz{}_dense_Re{}_Pm{}.pdf'.format(Nz, Reynolds, Pm)
else:
    fname = 'LOP_MA_scan_Nz{}_k{}_Re{}_Pm{}-test.pdf'.format(Nz, k, Reynolds, Pm)
ts = np.linspace(0.0, 80.0, 40)
Gs = np.zeros_like(ts)


def energy_norm_general(X1, X2, ma):
    u1 = X1['phi_z']
    w1 = -1.0j * kx * X1['phi']
    bx1 = X1['psi_z']
    bz1 = -1.0j * kx * X1['psi']
    u2 = X2['phi_z']
    w2 = -1.0j * kx * X2['phi']
    bx2 = X2['psi_z']
    bz2 = -1.0j * kx * X2['psi']

    field = (np.conj(u1)*u2 + np.conj(w1)*w2 + (1.0/ma**2.0)*(np.conj(bx1)*bx2 + np.conj(bz1)*bz2)).evaluate().integrate()
    return field['g'][0]


if Reynolds == np.inf:
    inviscid = True
else:
    inviscid = False
if mReynolds == np.inf:
    ideal = True
else:
    ideal = False
z_basis = de.Chebyshev('z', Nz, interval=(-0.5 * Lz, 0.5 * Lz))
domain = de.Domain([z_basis], grid_dtype=np.complex128)
probvars = ['phi', 'phi_z', 'psi', 'psi_z']
if not inviscid:
    probvars.append('uz')
    probvars.append('uzz')
else:  # really haven't tested this! May have neglected to add the right substitutions/equations
    probvars.append('zeta')
# problem = de.EVP(domain, variables=probvars, eigenvalue='omega')
problem = de.EVP(domain, variables=probvars, eigenvalue='gamma')
z = domain.grid(0)
a = 1.0

nccU = domain.new_field(name='U')
nccU['g'] = np.tanh(z / a)
problem.parameters['U'] = nccU

nccDU = domain.new_field(name='DU')
nccDU['g'] = 1 / (a * np.cosh(z / a) ** 2)
problem.parameters['DU'] = nccDU

problem.parameters['MA2'] = MA2 ** 2.0  # Alfven Mach number squared
problem.parameters['kx'] = kx
problem.parameters['Lz'] = Lz
if not inviscid:
    problem.parameters['Re'] = Reynolds
if not ideal:
    problem.parameters['Rm'] = mReynolds
problem.substitutions['dx(A)'] = "1.0j*kx*A"
# problem.substitutions['dt(A)'] = "-1.0j*omega*A"
problem.substitutions['dt(A)'] = "gamma*A"
problem.substitutions['u'] = 'phi_z'
problem.substitutions['w'] = "-dx(phi)"
problem.substitutions['Bz'] = "-dx(psi)"
problem.substitutions['Bx'] = "psi_z"
problem.substitutions['Jz'] = "dx(dx(psi)) + dz(Bx)"
if not inviscid:
    problem.substitutions['zeta'] = "dx(dx(phi)) + uz"
    problem.substitutions['zetaz'] = "dx(dx(u)) + uzz"
if not inviscid:
    problem.add_equation(
        "dt(zeta) - 1/Re*(dx(dx(zeta)) + dz(zetaz)) + U*dx(zeta) - dx(phi)*dz(DU) - 1/MA2*(dx(Jz)) = 0")
    problem.add_equation("uz - dz(u) = 0")
    problem.add_equation("uzz - dz(uz) = 0")
else:
    problem.add_equation("dt(zeta) + U*dx(zeta) - dx(phi)*dz(DU) - 1/MA2*(dx(Jz)) = 0")
    problem.add_equation("zeta - dx(dx(phi)) - dz(phi_z) = 0")
if not ideal:
    problem.add_equation("dt(psi) - 1/Rm*Jz + U*dx(psi) - dx(phi) = 0")
else:  # haven't tested this either
    problem.add_equation("dt(psi) + U*dx(psi) - dx(phi) = 0")
problem.add_equation("phi_z - dz(phi) = 0")
problem.add_equation("psi_z - dz(psi) = 0")

problem.add_bc("left(w) = 0")
problem.add_bc("right(w) = 0")
if not inviscid:
    problem.add_bc("left(dz(w)) = 0")  # no-slip, co-moving (with base flow U_0) boundaries
    problem.add_bc("right(dz(w)) = 0")  # (does free-slip make more sense?)
problem.add_bc("right(Bz) = 0")  # perfectly conducting boundaries
problem.add_bc("left(Bz) = 0")

logger.info("done with BCs")

EP = Eigenproblem(problem, grow_func=lambda x: x.imag, freq_func=lambda x: x.real, reject=False)  # TBD whether reject=False is the right move
pencil = 0

for ma_ind, ma in enumerate(MAs):
    if dense:
        EP.solve(sparse=False, pencil=pencil, parameters={'MA2': ma ** 2.0})
    else:
        EP.solve(sparse=True, N=k, pencil=pencil, parameters={'MA2': ma ** 2.0})
    pre_right = EP.solver.pencils[pencil].pre_right  # right-preconditioner for the matrix pencil, call it PR
    pre_right_LU = scipy.sparse.linalg.splu(pre_right.tocsc())  # LU factorization of the right-preconditioner
    if dense:
        finite_EVs = np.real(EP.solver.eigenvalues) < 10.0  # np.isfinite(EP.solver.eigenvalues)
        eigenvalues = EP.solver.eigenvalues[finite_EVs]
        eigenvectors = ((EP.solver.eigenvectors.T)[finite_EVs]).T  # pretty sure there's a better way to do this
    else:
        eigenvalues = EP.solver.eigenvalues
        eigenvectors = EP.solver.eigenvectors
    V = pre_right_LU.solve(eigenvectors)  # V satisfies PR @ V = solver.eigenvectors. I think this means V is the eigenvectors WITHOUT pre-conditioning
    Q, R = np.linalg.qr(V)  # Q is an orthonormal (in terms of the 2-norm!) basis in the subspace spanned by the eigenvectors in V, R is the change of basis tsfm
    M_E = EP.compute_mass_matrix(pre_right @ Q, lambda x1, x2: energy_norm_general(x1, x2, ma))  # This is M represented in the basis of the columns of Q
    Fdagger = np.linalg.cholesky(M_E)
    F = Fdagger.conj().T  # this gives ||F @ u||_2 = ||u||_E for u expressed in the (k-dimensional!) basis given by Q
    for tind, t in enumerate(ts):
        # K = R @ np.diag(np.exp(-1.0j*EP.solver.eigenvalues * t)) @ np.linalg.inv(R)
        K = R @ np.diag(np.exp(eigenvalues * t)) @ np.linalg.inv(R)
        s_vals = scipy.linalg.svdvals(F @ K @ np.linalg.inv(F))
        Gs[tind] = s_vals[0]
    Gs = Gs**2.0

    plt.semilogy(ts, Gs, c='C{}'.format(ma_ind), label=r'$M_A = {:4.2f}$'.format(ma))
    # plt.semilogy(ts, Gs[0]*np.exp(ts*2.0*np.max(np.imag(EP.evalues)))*Gs[-1]/np.exp(ts[-1]*2.0*np.max(np.imag(EP.evalues))), '--', c='C0')
    plt.semilogy(ts, Gs[0]*np.exp(ts*2.0*np.max(np.real(eigenvalues)))*Gs[-1]/np.exp(ts[-1]*2.0*np.max(np.real(eigenvalues))), '--', c='C{}'.format(ma_ind))
    print(ma_ind)
plt.xlabel(r'$T$')
plt.ylabel(r'$G(T)$')
plt.xlim((0, ts[-1]))
plt.legend()
plt.savefig(fname, bbox_inches='tight')
# plt.show()
