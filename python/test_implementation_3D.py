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

Nz = 256 # 512
Reynolds = 50.
Pm = 1.0
mReynolds = Pm * Reynolds
Lz = 10.0*np.pi
kx = 0.4
ky = 1e-8 # test 
MA = 0.5

k = 10
dense = True
reject_spurious = True

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

problem.parameters['MA2'] = MA ** 2.0  # Alfven Mach number squared
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
# problem.add_bc("right(Bz) = 0")
# problem.add_bc("right(dz(Bx)) = 0")
# problem.add_bc("right(dz(By)) = 0")
# problem.add_bc("left(Bz) = 0")
# problem.add_bc("left(dz(Bx)) = 0")
# problem.add_bc("left(dz(By)) = 0")

# problem.add_bc("right(Az) = 0")
problem.add_bc("right(Phi) = 0")
problem.add_bc("right(Ax) = 0")
problem.add_bc("right(Ay) = 0")
# problem.add_bc("left(Az) = 0")
problem.add_bc("left(Phi) = 0")
problem.add_bc("left(Ax) = 0")
problem.add_bc("left(Ay) = 0")

logger.info("done with BCs")

EP = Eigenproblem(problem, grow_func=lambda x: x.real, freq_func=lambda x: x.imag, reject=reject_spurious)  # TBD whether reject=False is the right move
pencil = 0

#target = -0.005-0.400j
target = 0
k = 4
EP.solve(sparse=True, N=k, pencil=pencil, target=target)
logger.info(f"sparse MA: {MA},  max eval = {np.max(np.real(EP.evalues))}")
indx = np.argmax(np.real(EP.evalues))
#fig = EP.plot_mode(0)
#fig.savefig("growth_mode_3D.png", dpi=300)
plot_vars = ['u','v','w','Bx','By']
with h5py.File(f"MA_{MA:.2f}_kx_{kx:.2f}_Re_{Reynolds:.0f}_Pm_{Pm:0.2f}_spectrum_3D.h5", "w") as df:
    df['spectrum'] = EP.evalues
    df['/scales/z'] = EP.grid()
    for m in range(k):
        emode = EP.eigenmode(m)
        for var in plot_vars:
            df[f'/mode_{m:d}/{var}'] = emode[var]['g'][:]
        # calculate Bz
        df[f'/mode_{m:d}/Bz'] = 1j*kx*emode['Ay']['g'][:] - 1j*ky*emode['Ax']['g'][:]
ax = EP.plot_spectrum()
ax.set_xlim(-6e-2, -2.5e-3)
ax.set_ylim(-0.45, 0.45)
ax.figure.savefig("spectrum_3D.png", dpi=300)
print(EP.evalues.real)
