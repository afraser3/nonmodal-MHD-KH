import h5py
from eigentools import Eigenproblem
from dedalus import public as de
import numpy as np
import scipy
from pathlib import Path
from matplotlib import pyplot as plt
import logging
plt.style.use('prl')

logger = logging.getLogger(__name__)

Nz = 256
Reynolds = 50.
Pm = 1.0
mReynolds = Pm * Reynolds
Lz = 10.0*np.pi
kx = 0.2
ky = 0.1
MA = 1e1

k = 100
dense = False
reject_spurious = True

z_basis = de.Chebyshev('z', Nz, interval=(-0.5 * Lz, 0.5 * Lz))
domain = de.Domain([z_basis], grid_dtype=np.complex128)
probvars = ['u', 'v', 'w', 'p', 'Bx', 'By', 'Bz', 'omega_x', 'omega_y', 'Bx_z', 'By_z', 'Bz_z']

problem = de.EVP(domain, variables=probvars, eigenvalue='gamma')
z = domain.grid(0)
a = 1.0

nccU = domain.new_field(name='U')
nccU['g'] = np.tanh(z / a)
# nccU['g'] = np.zeros_like(z)
problem.parameters['U'] = nccU

nccDU = domain.new_field(name='DU')
nccDU['g'] = 1 / (a * np.cosh(z / a) ** 2)
# nccDU['g'] = np.zeros_like(z)
problem.parameters['DU'] = nccDU

problem.parameters['MA2'] = MA ** 2  # Alfven Mach number squared
problem.parameters['kx'] = kx
problem.parameters['ky'] = ky
problem.parameters['Re'] = Reynolds
problem.parameters['Rm'] = mReynolds
problem.substitutions['dx(A)'] = "1j*kx*A"
problem.substitutions['dy(A)'] = "1j*ky*A"
problem.substitutions['dt(A)'] = "gamma*A"
problem.substitutions['Lap(A, Az)'] = '(dx(dx(A)) + dy(dy(A)) + dz(Az))'


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

problem.add_equation("dt(Bx) - dy(U*By - v) - dz(w - U*Bz) - Lap(Bx, Bx_z)/Rm = 0")
problem.add_equation("dt(By) - dx(U*By - v) - Lap(By, By_z)/Rm = 0")
problem.add_equation("dt(Bz) - dx(w - U*Bz) - Lap(Bz, Bz_z)/Rm = 0")
problem.add_equation("dz(Bx) - Bx_z = 0")
problem.add_equation("dz(By) - By_z = 0")
problem.add_equation("dz(Bz) - Bz_z = 0")

# no slip
problem.add_bc("left(u) = 0")
problem.add_bc("left(v) = 0")
problem.add_bc("left(w) = 0")
problem.add_bc("right(u) = 0")
problem.add_bc("right(v) = 0")
problem.add_bc("right(w) = 0")

# perfectly conducting boundaries
problem.add_bc("right(Bz) = 0")
# problem.add_bc("right(dz(Bx)) = 0")
problem.add_bc("right(Bx_z) = 0")
# problem.add_bc("right(dz(By)) = 0")
problem.add_bc("right(By_z) = 0")
problem.add_bc("left(Bz) = 0")
# problem.add_bc("left(dz(Bx)) = 0")
problem.add_bc("left(Bx_z) = 0")
# problem.add_bc("left(dz(By)) = 0")
problem.add_bc("left(By_z) = 0")

logger.info("done with BCs")


def div_B_int(X1):
    bx1 = X1['Bx']
    by1 = X1['By']
    bz_z1 = X1['Bz_z']
    field = (np.abs(1j*kx * bx1 + 1j*ky * by1 + bz_z1) ** 2).evaluate().integrate()
    return field['g'][0]


def energy_norm(X1):
    u1 = X1['u']
    v1 = X1['v']
    w1 = X1['w']
    bx1 = X1['Bx']
    by1 = X1['By']
    bz1 = X1['Bz']
    field = (np.abs(u1) ** 2 + np.abs(v1) ** 2 + np.abs(w1) ** 2 + (1.0 / MA**2.) * (np.abs(bx1) ** 2 + np.abs(by1) ** 2 + np.abs(bz1) ** 2)).evaluate().integrate()
    return field['g'][0]

EP = Eigenproblem(problem, grow_func=lambda x: x.real, freq_func=lambda x: x.imag, reject=reject_spurious)  # TBD whether reject=False is the right move
pencil = 0

target = 0
EP.solve(sparse=True, N=k, pencil=pencil, target=target)
logger.info(f"sparse MA: {MA},  max eval = {np.max(np.real(EP.evalues))},  num good evals = {len(EP.evalues)}")
indx = np.argmax(np.real(EP.evalues))
