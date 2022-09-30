import h5py
from eigentools import Eigenproblem
from dedalus import public as de
import numpy as np
import scipy
from pathlib import Path
from matplotlib import pyplot as plt
import logging
# plt.style.use('prl')

logger = logging.getLogger(__name__)

Nz = 256 # 512
Reynolds = 50.
Lz = 10.0*np.pi
kx = 0.25
ky = 0.25

k = 100
dense = False
reject_spurious = True

z_basis = de.Chebyshev('z', Nz, interval=(-0.5 * Lz, 0.5 * Lz))
domain = de.Domain([z_basis], grid_dtype=np.complex128)
probvars = ['u','v','w','p','omega_x','omega_y']

problem = de.EVP(domain, variables=probvars, eigenvalue='gamma')
z = domain.grid(0)
a = 1.0

nccU = domain.new_field(name='U')
nccU['g'] = np.tanh(z / a)
problem.parameters['U'] = nccU

nccDU = domain.new_field(name='DU')
nccDU['g'] = 1 / (a * np.cosh(z / a) ** 2)
problem.parameters['DU'] = nccDU

problem.parameters['kx'] = kx
problem.parameters['ky'] = ky
problem.parameters['Re'] = Reynolds
problem.substitutions['dx(A)'] = "1.0j*kx*A"
problem.substitutions['dy(A)'] = "1.0j*ky*A"
problem.substitutions['dt(A)'] = "gamma*A"

problem.substitutions['omega_z'] = "dx(v) - dy(u)"

problem.add_equation("dt(u) + U*dx(u) + w*DU + dx(p)          + (dy(omega_z) - dz(omega_y))/Re = 0")
problem.add_equation("dt(v) + U*dx(v)        + dy(p) + (dz(omega_x) - dx(omega_z))/Re = 0")
problem.add_equation("dt(w) + U*dx(w)        + dz(p) + (dx(omega_y) - dy(omega_x))/Re = 0")
problem.add_equation("dx(u) + dy(v) + dz(w) = 0")
problem.add_equation("dy(w) - dz(v) - omega_x = 0")
problem.add_equation("dz(u) - dx(w) - omega_y = 0")

# no slip
problem.add_bc("left(u) = 0")
problem.add_bc("left(v) = 0")
problem.add_bc("left(w) = 0")
problem.add_bc("right(u) = 0")
problem.add_bc("right(v) = 0")
problem.add_bc("right(w) = 0")

logger.info("done with BCs")

EP = Eigenproblem(problem, grow_func=lambda x: x.real, freq_func=lambda x: x.imag, reject=reject_spurious)  # TBD whether reject=False is the right move
pencil = 0

#target = -0.005-0.400j
target = 0
EP.solve(sparse=True, N=k, pencil=pencil, target=target)
# logger.info(f"sparse MA: {MA},  max eval = {np.max(np.real(EP.evalues))}")
# indx = np.argmax(np.real(EP.evalues))
#fig = EP.plot_mode(0)
#fig.savefig("growth_mode_3D.png", dpi=300)
# plot_vars = ['u','v','w','Bx','By']
# with h5py.File(f"MA_{MA:.2f}_kx_{kx:.2f}_Re_{Reynolds:.0f}_Pm_{Pm:0.2f}_spectrum_3D.h5", "w") as df:
    # df['spectrum'] = EP.evalues
    # df['/scales/z'] = EP.grid()
    # for m in range(k):
        # emode = EP.eigenmode(m)
        # for var in plot_vars:
            # df[f'/mode_{m:d}/{var}'] = emode[var]['g'][:]
        # calculate Bz
        # df[f'/mode_{m:d}/Bz'] = 1j*kx*emode['Ay']['g'][:] - 1j*ky*emode['Ax']['g'][:]
# ax = EP.plot_spectrum()
# ax.set_xlim(-6e-2, -2.5e-3)
# ax.set_ylim(-0.45, 0.45)
# ax.figure.savefig("spectrum_3D.png", dpi=300)
# print(EP.evalues.real)
