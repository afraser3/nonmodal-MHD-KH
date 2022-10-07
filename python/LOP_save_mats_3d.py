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
ky = 0.025
# MA = 10.0  # Alfven Mach number. For this value, the system is unstable.
MA = 1.0
MA2 = MA**2.0
k = 300
dense = False
reject_spurious = True
if dense:
    savedir = Path("LOPs_3D/individual/Re{}/MA{}/Pm{}/kx{}/ky{}".format(Reynolds, MA, Pm, kx, ky))
    fname = str(savedir / 'Nz{}_dense.h5'.format(Nz))
    # fname = "LOPs_3D/individual/Re{}/MA{}/Pm{}/2d/Nz{}_dense.h5".format(Reynolds, MA, Pm, Nz)
else:
    savedir = Path("LOPs_3D/individual/Re{}/MA{}/Pm{}/kx{}/ky{}".format(Reynolds, MA, Pm, kx, ky))
    fname = str(savedir / 'Nz{}_k{}.h5'.format(Nz, k))
    # fname = "LOPs_3D/individual/Re{}/MA{}/Pm{}/Nz{}_k{}.h5".format(Reynolds, MA, Pm, Nz, k)
if not savedir.exists():
    savedir.mkdir(parents=True)
# ts = np.linspace(0.0, 50.0, 51)
ts = np.linspace(0.0, 200.0, 201)
# ts = np.linspace(0.0, 1000.0, 501)


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

EP = Eigenproblem(problem, grow_func=lambda x: x.real, freq_func=lambda x: x.imag, reject=reject_spurious)
pencil = 0

Gs, eigenvalues, eigenvectors, F, Q, R, pre_right = LOP(EP, ts, {'MA2': MA2}, lambda x1, x2: energy_norm_general(x1, x2, MA), pencil=pencil, dense=dense, k=k, return_setup_mats=True)
Finv = np.linalg.inv(F)

with h5py.File(fname, mode="w-") as file:
    setup_mats_grp = file.create_group('setup_mats')
    evalues_data = setup_mats_grp.create_dataset('eigenvalues', data=eigenvalues)
    evecs_data = setup_mats_grp.create_dataset('eigenvectors', data=eigenvectors)
    R_data = setup_mats_grp.create_dataset('R', data=R)
    Q_data = setup_mats_grp.create_dataset('Q', data=Q)
    F_data = setup_mats_grp.create_dataset('F', data=F)
    pre_right_data = setup_mats_grp.create_dataset('pre_right', data=pre_right.toarray())

    scales_grp = file.create_group('scales')  # in the spirit of "scales" in d2 h5 files, but this is a bad name for this use case
    ts_data = scales_grp.create_dataset('ts', data=ts)
    kx_data = scales_grp.create_dataset('kx', data=kx)
    ky_data = scales_grp.create_dataset('ky', data=ky)
    Nz_data = scales_grp.create_dataset('Nz', data=Nz)
    Lz_data = scales_grp.create_dataset('Lz', data=Lz)
    zs_data = scales_grp.create_dataset('z', data=z)
    print("starting time loop")

    energies_grp = file.create_group('energies')  # maybe I'll put ME or KE or something here, not just Gs?
    Gs_data = energies_grp.create_dataset('Gs', data=Gs)

    for tind, t in enumerate(ts):
        # if False:
        t_grp = file.create_group('ti_{}'.format(tind))
        K = R @ np.diag(np.exp(eigenvalues * t)) @ np.linalg.inv(R)
        propagator = F @ K @ Finv
        U, s, Vh = scipy.linalg.svd(propagator)
        # prop_data = t_grp.create_dataset('propagator', data=propagator)
        # s_data = t_grp.create_dataset('s', data=s)
        # Vh_data = t_grp.create_dataset('Vh', data=Vh)
        v1 = Vh.conj()[0]  # the optimal perturbation in the F'ed and Q'ed basis
        v1F = Finv @ v1  # change of basis away from the F'ed basis
        # now I want to say that R is the change-of-basis matrix between the Q and eigenvector bases, and so I should
        # be doing R^-1, but that doesn't give the right dimensions.
        # I'm pretty sure that v1F is v1 in the basis of the columns of Q, so I should be able to do Q @ v1F to
        # get the perturbation in terms of Chebyshev expansion coefficients
        LOP = Q @ v1F
        LOP = pre_right @ LOP
        LOP_phase = LOP[0] / np.abs(LOP[0])
        LOP = LOP / LOP_phase
        LOP_data = t_grp.create_dataset('LOP', data=LOP)
        v1_data = t_grp.create_dataset('v1', data=v1)
print("done")
