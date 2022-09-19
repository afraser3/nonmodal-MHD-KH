import h5py
from eigentools import Eigenproblem
from dedalus import public as de
import numpy as np
import scipy
from pathlib import Path
from matplotlib import pyplot as plt
import logging
logger = logging.getLogger(__name__)
plt.style.use('prl')

Nz = 256 #512
Reynolds = 50.0
Pm = 1.0
mReynolds = Pm * Reynolds
Lz = 10.0*np.pi
kx = 0.4
# MA = 10.0  # Alfven Mach number. For this value, the system is unstable.
MA = 0.5  # For this value, the system is stable.
MA2 = MA**2.0
k = 2
dense = True
reject_spurious = True

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

problem.parameters['MA2'] = MA2 # Alfven Mach number squared
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

EP = Eigenproblem(problem, grow_func=lambda x: x.real, freq_func=lambda x: x.imag, reject=reject_spurious)  # TBD whether reject=False is the right move
pencil = 0

sp_growth, sp_indx, sp_freq = EP.growth_rate(sparse=True)
logger.info(f"sparse MA: {MA},  max eval = {EP.evalues[sp_indx]}")
fig = EP.plot_mode(sp_indx)
fig.savefig("growth_mode_2D.png", dpi=300)

ax = EP.plot_spectrum()
ax.set_xlim(-6e-2, -2.5e-3)
ax.set_ylim(-0.45, 0.45)
ax.figure.savefig("spectrum_2D.png", dpi=300)
print(EP.evalues.real)
with h5py.File(f"MA_{MA:.2f}_kx_{kx:.2f}_Re_{Reynolds:.0f}_Pm_{Pm:0.2f}_spectrum_2D.h5", "w") as df:
    df['spectrum'] = EP.evalues
    spec_size = df['spectrum'].shape[0]
    df['/scales/z'] = EP.grid()
    for m in range(spec_size):
        emode = EP.eigenmode(m)
        df[f'/mode_{m:d}/u'] = emode['phi_z']['g'][:]
        df[f'/mode_{m:d}/w'] = -1j*kx*emode['phi']['g'][:]
        df[f'/mode_{m:d}/Bx'] = emode['psi_z']['g'][:]
        df[f'/mode_{m:d}/Bz'] = -1j*kx*emode['psi']['g'][:]
