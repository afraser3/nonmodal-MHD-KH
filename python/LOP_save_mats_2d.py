from linear_optimals import LOP_setup_mats, LOP, F_from_norm

from eigentools import Eigenproblem
from dedalus import public as de
import numpy as np
import scipy
from pathlib import Path
import h5py
import logging
logger = logging.getLogger(__name__)

Nz = 256  # 128 is definitely not enough for these parameters
Reynolds = 50.0
Pm = 1.0
mReynolds = Pm * Reynolds
Lz = 10.0*np.pi
kx = 0.2
# MA = 10.0  # Alfven Mach number. For this value, the system is unstable.
MA = 0.5
MA2 = MA**2.0
k = 100
dense = False
reject_spurious = True
if dense:  # TODO: the name should really include kx
    savedir = Path("LOPs/Re{}/MA{}/Pm{}/2d".format(Reynolds, MA, Pm))
    fname = "LOPs/Re{}/MA{}/Pm{}/2d/Nz{}_dense.h5".format(Reynolds, MA, Pm, Nz)
else:
    savedir = Path("LOPs/Re{}/MA{}/Pm{}/2d".format(Reynolds, MA, Pm))
    fname = "LOPs/Re{}/MA{}/Pm{}/2d/Nz{}_k{}.h5".format(Reynolds, MA, Pm, Nz, k)
if not savedir.exists():
    savedir.mkdir(parents=True)
# ts = np.linspace(0.0, 50.0, 51)
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


def KE_norm_general(X1, X2, ma):  # ma should not be an argument here, but I'm not sure exactly how to generalize the relevant part of linear_optimals.F_from_norm
    u1 = X1['phi_z']
    w1 = -1.0j * kx * X1['phi']
    u2 = X2['phi_z']
    w2 = -1.0j * kx * X2['phi']

    field = (np.conj(u1)*u2 + np.conj(w1)*w2).evaluate().integrate()
    return field['g'][0]


def ME_norm_general(X1, X2, ma):
    bx1 = X1['psi_z']
    bz1 = -1.0j * kx * X1['psi']
    bx2 = X2['psi_z']
    bz2 = -1.0j * kx * X2['psi']

    field = ((1.0/ma**2.0)*(np.conj(bx1)*bx2 + np.conj(bz1)*bz2)).evaluate().integrate()
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

EP = Eigenproblem(problem, grow_func=lambda x: x.imag, freq_func=lambda x: x.real, reject=reject_spurious)
pencil = 0

eigenvalues, eigenvectors, F, Q, R, pre_right = LOP_setup_mats(EP, {'MA2': MA2}, energy_norm_general, MA, pencil=pencil, dense=dense)
Gs, eigenvalues, eigenvectors = LOP(EP, ts, {'MA2': MA2}, energy_norm_general, MA, pencil=pencil, dense=dense)
Finv = np.linalg.inv(F)
# F_KE = F_from_norm(EP, pre_right, Q, KE_norm_general, MA)
# F_KE_inv = np.linalg.inv(F_KE)
# F_ME = F_from_norm(EP, pre_right, Q, ME_norm_general, MA)
# F_ME_inv = np.linalg.inv(F_ME)

with h5py.File(fname, mode="w-") as file:
    setup_mats_grp = file.create_group('setup_mats')
    evalues_data = setup_mats_grp.create_dataset('eigenvalues', data=eigenvalues)
    evecs_data = setup_mats_grp.create_dataset('eigenvectors', data=eigenvectors)
    R_data = setup_mats_grp.create_dataset('R', data=R)
    Q_data = setup_mats_grp.create_dataset('Q', data=Q)
    F_data = setup_mats_grp.create_dataset('F', data=F)
    # F_KE_data = setup_mats_grp.create_dataset('F_KE', data=F_KE)
    # F_ME_data = setup_mats_grp.create_dataset('F_ME', data=F_ME)
    pre_right_data = setup_mats_grp.create_dataset('pre_right', data=pre_right.toarray())

    scales_grp = file.create_group('scales')  # in the spirit of "scales" in d2 h5 files, but this is a bad name for this use case
    ts_data = scales_grp.create_dataset('ts', data=ts)
    kx_data = scales_grp.create_dataset('kx', data=kx)
    Nz_data = scales_grp.create_dataset('Nz', data=Nz)
    Lz_data = scales_grp.create_dataset('Lz', data=Lz)
    zs_data = scales_grp.create_dataset('z', data=z)
    print("starting time loop")

    energies_grp = file.create_group('energies')  # maybe I'll put ME or KE or something here, not just Gs?
    Gs_data = energies_grp.create_dataset('Gs', data=Gs)

    for tind, t in enumerate(ts):
        t_grp = file.create_group('ti_{}'.format(tind))
        # K = R @ np.diag(np.exp(-1.0j*EP.solver.eigenvalues * t)) @ np.linalg.inv(R)
        K = R @ np.diag(np.exp(eigenvalues * t)) @ np.linalg.inv(R)
        propagator = F @ K @ Finv
        U, s, Vh = scipy.linalg.svd(propagator)
        # s_vals = scipy.linalg.svdvals(F @ K @ np.linalg.inv(F))
        # Gs[tind] = s_vals[0]**2.0
        # K_data = t_grp.create_dataset('K', data=K)
        prop_data = t_grp.create_dataset('propagator', data=propagator)
        # U_data = t_grp.create_dataset('U', data=U)
        s_data = t_grp.create_dataset('s', data=s)
        Vh_data = t_grp.create_dataset('Vh', data=Vh)
print("done")
