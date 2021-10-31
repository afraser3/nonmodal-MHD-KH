from eigentools import Eigenproblem
from dedalus import public as de
import numpy as np
import scipy
from matplotlib import pyplot as plt
import logging
logger = logging.getLogger(__name__)

Nz = 256
Reynolds = 50.0
mReynolds = 50.0
Lz = 10.0*np.pi
kx = 0.2
MA = 0.1  # 25.0
MA2 = MA**2.0
k = 100
ts = np.linspace(0.0, 30.0, 20)
Gs = np.zeros_like(ts)


def energy_norm_general(X1, X2, MA):  # because the norm depends on MA, need to include MA as an argument
    u1 = X1['phi_z']
    w1 = -1.0j * kx * X1['phi']
    bx1 = X1['psi_z']
    bz1 = -1.0j * kx * X1['psi']
    u2 = X2['phi_z']
    w2 = -1.0j * kx * X2['phi']
    bx2 = X2['psi_z']
    bz2 = -1.0j * kx * X2['psi']

    field = (np.conj(u1)*u2 + np.conj(w1)*w2 + (1.0/MA**2.0)*(np.conj(bx1)*bx2 + np.conj(bz1)*bz2)).evaluate().integrate()
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
problem = de.EVP(domain, variables=probvars, eigenvalue='omega')
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
problem.substitutions['dt(A)'] = "-1.0j*omega*A"
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

EP = Eigenproblem(problem, grow_func=lambda x: x.imag, freq_func=lambda x: x.real)

#k, (real_points, imag_points), mu=mu, inner_product=lambda x1, x2: energy_norm_general(x1, x2, ma)
# inner_product = lambda x1, x2: energy_norm_general(x1, x2, MA)
pencil = 0
EP.solve(sparse=True, N=k, pencil=pencil)
pre_right = EP.solver.pencils[pencil].pre_right  # right-preconditioner for the matrix pencil, call it PR
pre_right_LU = scipy.sparse.linalg.splu(pre_right.tocsc())  # LU factorization of the right-preconditioner
V = pre_right_LU.solve(EP.solver.eigenvectors)  # V satisfies PR @ V = solver.eigenvectors. I think this means V is the eigenvectors WITHOUT pre-conditioning
Q, R = np.linalg.qr(V)
E = -(EP.solver.pencils[pencil].M_exp.toarray())
E_inv = scipy.linalg.pinv(E)  # this is probably bad
A = (EP.solver.pencils[pencil].L_exp.toarray())
# E = -(EP.solver.pencils[pencil].M_exp)
# E_inv = scipy.sparse.linalg.inv(E)
# A = (EP.solver.pencils[pencil].L_exp)
E_inv_A = E_inv @ A
E_inv_A_Qbasis = np.conj(Q.T) @ E_inv_A @ Q
M = EP.compute_mass_matrix(pre_right@Q, lambda x1, x2: energy_norm_general(x1, x2, MA))
F = scipy.linalg.cholesky(M)
F_inv = scipy.linalg.inv(F)
for tind, t in enumerate(ts):
    # K = scipy.linalg.expm(E_inv @ A * t)
    K = scipy.linalg.expm(E_inv_A_Qbasis * t)
    s_vals = scipy.linalg.svdvals(F@K@F_inv)
    Gs[tind] = s_vals[0]
plt.semilogy(ts, Gs)
plt.semilogy(ts, Gs[0]*np.exp(ts*np.max(np.imag(EP.evalues)))*Gs[-1]/np.exp(ts[-1]*np.max(np.imag(EP.evalues))), '--', c='C0')
plt.xlabel(r'$T$')
plt.ylabel(r'$G(T)$')
plt.show()
