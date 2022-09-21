from linear_optimals import LOP

from eigentools import Eigenproblem
from dedalus import public as de
import numpy as np
import scipy
from pathlib import Path
from matplotlib import pyplot as plt
import logging
logger = logging.getLogger(__name__)
# plt.style.use('prl')

Nz = 256
Reynolds = 50.0
Pm = 1.0
mReynolds = Pm * Reynolds
Lz = 10.0*np.pi
kx = 0.2
# MA = 10.0  # Alfven Mach number. For this value, the system is unstable.
MA = 1.0  # For this value, the system is stable.
MAs = np.array([0.5, 0.75, 1.0])
MA2 = MA**2.0
k = 100
dense = False
reject_spurious = True
plotdir = Path("LOPs/Re{}".format(Reynolds))
if not plotdir.exists():
    plotdir.mkdir(parents=True)
if dense:
    fname = 'LOPs/Re{}/LOP_MA_scan_kx{}_Nz{}_dense_Pm{}_2d.pdf'.format(Reynolds, kx, Nz, Pm)
else:
    fname = 'LOPs/Re{}/LOP_MA_scan_kx{}_Nz{}_k{}_Pm{}_2d.pdf'.format(Reynolds, kx, Nz, k, Pm)
ts = np.linspace(0.0, 80.0, 200)


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
problem = de.EVP(domain, variables=probvars, eigenvalue='gamma')
z = domain.grid(0)
a = 1.0

nccU = domain.new_field(name='U')
nccU['g'] = np.tanh(z / a)
problem.parameters['U'] = nccU

nccDU = domain.new_field(name='DU')
nccDU['g'] = 1 / (a * np.cosh(z / a) ** 2)
problem.parameters['DU'] = nccDU

problem.parameters['MA2'] = MA2  # Alfven Mach number squared
problem.parameters['kx'] = kx
problem.parameters['Lz'] = Lz
if not inviscid:
    problem.parameters['Re'] = Reynolds
if not ideal:
    problem.parameters['Rm'] = mReynolds
problem.substitutions['dx(A)'] = "1.0j*kx*A"
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

EP = Eigenproblem(problem, grow_func=lambda x: x.real, freq_func=lambda x: x.imag, reject=reject_spurious)  # TBD whether reject=False is the right move
pencil = 0

for ma_ind, ma in enumerate(MAs):
    Gs, eigenvalues, eigenvectors = LOP(EP, ts, {'MA2': ma ** 2.0}, energy_norm_general, ma, pencil=pencil, dense=dense)

    plt.semilogy(ts, Gs, c='C{}'.format(ma_ind), label=r'$M_A = {:4.2f}$'.format(ma))
    plt.semilogy(ts, Gs[0]*np.exp(ts*2.0*np.max(np.real(eigenvalues)))*Gs[-1]/np.exp(ts[-1]*2.0*np.max(np.real(eigenvalues))), '--', c='C{}'.format(ma_ind))
    logger.info(f"MA: {ma},  max eval = {np.max(np.real(eigenvalues))}")
plt.xlabel(r'$T$')
plt.ylabel(r'$G(T)$')
plt.xlim((0, ts[-1]))
plt.legend()
plt.savefig(fname, bbox_inches='tight')
