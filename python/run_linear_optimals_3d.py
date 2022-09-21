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
Reynolds = 50.
Pm = 1.0
mReynolds = Pm * Reynolds
Lz = 10.0*np.pi
kx = 0.2
ky = 1e-8  # worked for ky=0.1, but 0.2 gives no good eigenvalues, strangely. At Nz=1024 it gave 1 good eigenvalue.
# MA = 10.0  # Alfven Mach number. For this value, the system is unstable.
MA = 1.0  # For this value, the system is stable.
MAs = np.array([0.5, 0.75, 1.0])
MA2 = MA**2.0
k = 100
dense = False
reject_spurious = True
plotdir = Path("LOPs_3D/Re{}".format(Reynolds))
if not plotdir.exists():
    plotdir.mkdir(parents=True)
if dense:
    fname = 'LOPs_3D/Re{}/LOP_MA_scan_kx{}_ky{}_Nz{}_dense_Pm{}_3d.pdf'.format(Reynolds, kx, ky, Nz, Pm)
else:
    fname = 'LOPs_3D/Re{}/LOP_MA_scan_kx{}_ky{}_Nz{}_k{}_Pm{}_3d.pdf'.format(Reynolds, kx, ky, Nz, k, Pm)
ts = np.linspace(0.0, 80.0, 200)


def energy_norm_general(X1, X2, ma):
    u1 = X1['u']
    v1 = X1['v']
    w1 = X1['w']
    bx1 = X1['Bx']
    by1 = X1['By']
    bz1 = 1j * kx * X1['Ay'] - 1j * ky * X1['Ax']

    u2  = X2['u']
    v2  = X2['v']
    w2  = X2['w']
    bx2 = X2['Bx']
    by2 = X2['By']
    bz2 = 1j * kx * X2['Ay'] - 1j * ky * X2['Ax']

    field = (np.conj(u1)*u2 + np.conj(v1)*v2 + np.conj(w1)*w2 + (1.0/ma**2.0)*(np.conj(bx1)*bx2 + np.conj(by1)*by2 + np.conj(bz1)*bz2)).evaluate().integrate()
    return field['g'][0]

z_basis = de.Chebyshev('z', Nz, interval=(-0.5 * Lz, 0.5 * Lz))
domain = de.Domain([z_basis], grid_dtype=np.complex128)
probvars = ['u','v','w','p','Ax','Ay','Az','Phi','omega_x','omega_y','Bx','By']

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
problem.parameters['ky'] = ky
problem.parameters['Re'] = Reynolds
problem.parameters['Rm'] = mReynolds
problem.substitutions['dx(A)'] = "1.0j*kx*A"
problem.substitutions['dy(A)'] = "1.0j*ky*A"
problem.substitutions['dt(A)'] = "gamma*A"

problem.substitutions['Bz'] = "dx(Ay) - dy(Ax)"
problem.substitutions['omega_z'] = "dx(v) - dy(u)"

# Coulumb gauge 
problem.substitutions['Jx'] = "dy(Bz) - dz(By)"
problem.substitutions['Jy'] = "dz(Bx) - dx(Bz)"
problem.substitutions['Jz'] = "dx(By) - dy(Bx)"

problem.add_equation("dt(u) + U*dx(u) + w*DU + dx(p)          + (dy(omega_z) - dz(omega_y))/Re = 0")
problem.add_equation("dt(v) + U*dx(v)        + dy(p) - Jz/MA2 + (dz(omega_x) - dx(omega_z))/Re = 0")
problem.add_equation("dt(w) + U*dx(w)        + dz(p) + Jy/MA2 + (dx(omega_y) - dy(omega_x))/Re = 0")
problem.add_equation("dx(u) + dy(v) + dz(w) = 0")
problem.add_equation("dy(w) - dz(v) - omega_x = 0")
problem.add_equation("dz(u) - dx(w) - omega_y = 0")

problem.add_equation("dt(Ax) + dx(Phi)             + Jx/Rm = 0")
problem.add_equation("dt(Ay) + dy(Phi) - w + U*Bz  + Jy/Rm = 0")
problem.add_equation("dt(Az) + dz(Phi) + v - U*By  + Jz/Rm = 0")
problem.add_equation("dx(Ax) + dy(Ay) + dz(Az) = 0")
problem.add_equation("dy(Az) - dz(Ay) - Bx  = 0")
problem.add_equation("dz(Ax) - dx(Az) - By  = 0")

# no slip
problem.add_bc("left(u) = 0")
problem.add_bc("left(v) = 0")
problem.add_bc("left(w) = 0")
problem.add_bc("right(u) = 0")
problem.add_bc("right(v) = 0")
problem.add_bc("right(w) = 0")

# perfectly conducting boundaries
problem.add_bc("right(Bz) = 0")
problem.add_bc("right(dz(Bx)) = 0")
problem.add_bc("right(dz(By)) = 0")
problem.add_bc("left(Bz) = 0")
problem.add_bc("left(dz(Bx)) = 0")
problem.add_bc("left(dz(By)) = 0")

logger.info("done with BCs")

EP = Eigenproblem(problem, grow_func=lambda x: x.real, freq_func=lambda x: x.imag, reject=reject_spurious)  # TBD whether reject=False is the right move
pencil = 0

for ma_ind, ma in enumerate(MAs):
    Gs, eigenvalues, eigenvectors = LOP(EP, ts, {'MA2': ma ** 2.0}, energy_norm_general, ma, pencil=pencil, dense=dense, k=k)

    plt.semilogy(ts, Gs, c='C{}'.format(ma_ind), label=r'$M_A = {:4.2f}$'.format(ma))
    plt.semilogy(ts, Gs[0]*np.exp(ts*2.0*np.max(np.real(eigenvalues)))*Gs[-1]/np.exp(ts[-1]*2.0*np.max(np.real(eigenvalues))), '--', c='C{}'.format(ma_ind))
    logger.info(f"MA: {ma},  max eval = {np.max(np.real(eigenvalues))}")
plt.xlabel(r'$T$')
plt.ylabel(r'$G(T)$')
plt.xlim((0, ts[-1]))
plt.legend()
plt.savefig(fname, bbox_inches='tight')
