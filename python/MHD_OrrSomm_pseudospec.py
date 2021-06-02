"""
Copying this from an old script I had lying around from grad school.
I did very little testing of the viscous/resistive system, so it needs to be checked.
(I have verified that it returns a reasonable eigenvalue for the current set of parameters.)
It could definitely also use cleaning up.
"""
from eigentools import Eigenproblem
from dedalus import public as de
import numpy as np
from matplotlib import pyplot as plt

Reynolds = 10.0
mReynolds = 1.0
Pm = mReynolds/Reynolds

Nz = 512
Lz = 10.0*np.pi
kx = 0.01
MA = 0.99


def energy_norm(X1, X2):
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
else:
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

problem.parameters['MA2'] = MA ** 2.0  # Alfven Mach number squared
problem.parameters['kx'] = kx
problem.parameters['Lz'] = Lz
if not inviscid:
    problem.parameters['Re'] = Reynolds
if not ideal:
    problem.parameters['Rm'] = mReynolds
problem.substitutions['dx(A)'] = "1.0j*kx*A"
# problem.substitutions['dt(A)'] = "1.0j*omega*A"
# problem.substitutions['dt(A)'] = "sigma*A"  # was having trouble with EP.calc_ps using omega
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
else:
    problem.add_equation("dt(psi) + U*dx(psi) - dx(phi) = 0")
problem.add_equation("phi_z - dz(phi) = 0")
problem.add_equation("psi_z - dz(psi) = 0")

problem.add_bc("left(w) = 0")
problem.add_bc("right(w) = 0")
if not inviscid:
    problem.add_bc("left(dz(w)) = 0")  # no-slip, co-moving (with base flow U_0) boundaries
    problem.add_bc("right(dz(w)) = 0")
problem.add_bc("right(Bz) = 0")  # perfectly conducting boundaries
problem.add_bc("left(Bz) = 0")

print("done with BCs")

EP = Eigenproblem(problem, grow_func=lambda x: x.imag, freq_func=lambda x: x.real)
# EP.solve(sparse=False)
# print("done with EP.solve()")

k = 100

psize = 100
real_points = np.linspace(-0.5, 0.5, psize)
imag_points = np.linspace(-0.5, 0.5, psize)
EP.calc_ps(k, (real_points, imag_points), inner_product=energy_norm)

plt.plot(np.real(EP.evalues_low), np.imag(EP.evalues_low), 'x', label=r'evalues_low')
plt.plot(np.real(EP.evalues_high), np.imag(EP.evalues_high), '+', label=r'evalues_high')
if len(EP.evalues) > 0:
    plt.plot(np.real(EP.evalues), np.imag(EP.evalues), '.', label=r'evalues')
plt.contour(EP.ps_real, EP.ps_imag, np.log10(EP.pseudospectrum), levels=np.arange(-8, 0))
plt.colorbar(label=r'$\log_{10} (\epsilon)$')
plt.legend()
plt.ylim((-0.5, 0.5))
plt.xlim((-0.5, 0.5))
plt.axhline(0, color='k', alpha=0.2)
plt.xlabel('real (oscillating) frequency')
plt.ylabel('growth rate')
plt.title(r'(kx, MA, Re, Pm) = ({}, {}, {}, {})'.format(kx, MA, Reynolds, Pm))
plt.savefig('MHD_KH_pseudospectra_kx{}_MA{}_Re{}_Pm{}_Nz{}.pdf'.format(kx, MA, Reynolds, Pm, Nz))
plt.show()
