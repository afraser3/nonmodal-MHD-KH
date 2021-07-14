"""
Calculates the epsilon-pseudospectrum for a shear layer in 2D incompressible MHD.
Saves results to an hdf5 file.
Parameters are hard-coded.
TODO: learn how to use ConfigParser so I don't have to mess with git committing hard-coded parameters.
TODO: MPI-parallelize
"""
from eigentools import Eigenproblem
from dedalus import public as de
import numpy as np
import os
import h5py

Reynolds = 1000.0
mReynolds = 1000.0
kx = 0.2
MA = 1.2

Pm = mReynolds/Reynolds
Nz = 512
Lz = 10.0*np.pi

k = 200
psize = 100
freq_axis_bounds = [-1.0 * kx / 0.4, 1.0 * kx / 0.4]  # the reasonable/helpful bounds seem to scale roughly with kx
growth_axis_bounds = [-1.0 * kx / 0.4, 0.2]  # 0.2 is sometimes too small to see the full extent, but usually fine
real_points = np.linspace(freq_axis_bounds[0], freq_axis_bounds[1], psize)
imag_points = np.linspace(growth_axis_bounds[0], growth_axis_bounds[1], psize)

filepath = 'saved_spectra/Pm{}/Re{}/MA{}/kx{}/'.format(Pm, Reynolds, MA, kx)
filename = 'pseudospec_Nz{}_Lz{}pi_k{}.h5'.format(Nz, Lz/np.pi, k)

try:
    os.makedirs(filepath)
except FileExistsError:
    pass
if filename in os.listdir(filepath):
    raise FileExistsError

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

EP.calc_ps(k, (real_points, imag_points), inner_product=energy_norm)

with h5py.File(filepath+filename, 'w-') as file:
    evalues_grp = file.create_group('evalues')
    evecs_grp = file.create_group('evecs')
    pseudospec_grp = file.create_group('pseudospec')

    evalues = evalues_grp.create_dataset('evalues', data=EP.evalues)
    evalues_low = evalues_grp.create_dataset('evalues_low', data=EP.evalues_low)
    evalues_high = evalues_grp.create_dataset('evalues_high', data=EP.evalues_high)

    evectors = evecs_grp.create_dataset('evectors', data=EP.solver.eigenvectors)

    ps_real = pseudospec_grp.create_dataset('ps_real', data=EP.ps_real)
    ps_imag = pseudospec_grp.create_dataset('ps_imag', data=EP.ps_imag)
    pseudospectrum = pseudospec_grp.create_dataset('pseudospectrum', data=EP.pseudospectrum)
    # I think ps_real, ps_imag are the same thing as real_points and imag_points, and thus it would be redundant
    # to save those as well. Right?
