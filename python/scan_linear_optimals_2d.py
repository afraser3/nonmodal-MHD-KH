"""
Does all the stuff in LOP_save_mats but as a kx, ky scan.

To run in parallel with 2 processes, e.g., do:
    $ mpiexec -n 2 python3 scan_linear_optimals_2d.py LOPs_2D/config_files/run1.cfg
"""
from linear_optimals import LOP

from eigentools import Eigenproblem
from dedalus import public as de
import numpy as np
import scipy
import sys
from pathlib import Path
import h5py
from matplotlib import pyplot as plt
import logging
from configparser import ConfigParser
from mpi4py import MPI
CW = MPI.COMM_WORLD
logger = logging.getLogger(__name__)
# plt.style.use('prl')

# Parse the .cfg filename passed to script
config_file = Path(sys.argv[-1])
logger.info("Running with config file {}".format(str(config_file)))

# Parse settings from .cfg file and save to variables
runconfig = ConfigParser()
runconfig.read(str(config_file))
# data dir to save to: if "LOPs_2D/config_files/run1.cfg" is input, then "LOPs_2D/run1/" becomes the output directory
basedir = Path("LOPs_2D") / config_file.stem
if CW.rank == 0:  # only do this for the 0th MPI process
    if not basedir.exists():
        basedir.mkdir(parents=True)

# Parameters
params = runconfig['parameters']
MA = params.getfloat('MA')
MA2 = MA**2.0
Reynolds = params.getfloat('Reynolds')
Pm = params.getfloat('Pm')
mReynolds = Pm*Reynolds
kx_start = params.getfloat('kx_start')
kx_stop = params.getfloat('kx_stop')
kx_num = params.getint('kx_num')
kx_endpoint = params.getboolean('kx_endpoint')

t_start = params.getfloat('t_start')
t_stop = params.getfloat('t_stop')
t_num = params.getint('t_num')
t_endpoint = params.getboolean('t_endpoint')
ts = np.linspace(t_start, t_stop, t_num, endpoint=t_endpoint)

k = params.getint('k')
dense = params.getboolean('dense')
Nz = params.getint('Nz')
Lz = params.getfloat('Lz_factor')*np.pi
reject_spurious = True
make_plots = params.getboolean('make_plots')

kx_global = np.linspace(kx_start, kx_stop, kx_num, endpoint=kx_endpoint)  # the full set of kxs to scan over
kx_local = kx_global[CW.rank::CW.size]  # the kxs that this processor will scan over
logger.debug('This MPI process will scan over these kxs: {}'.format(kx_local))

kxdirs_global = [basedir / 'kxi_{:04d}/'.format(kx_i) for kx_i, kx in enumerate(kx_global)]
kxdirs_local = kxdirs_global[CW.rank::CW.size]
if CW.rank == 0:
    for kxdir in kxdirs_global:
        if not kxdir.exists():
            kxdir.mkdir(parents=True)

kx0 = kx_global[0]  # just a placeholder MA for setting Dedalus EVP up
if CW.rank == 0:
    logger.info('Scanning over these kxs: {}'.format(kx_global))

# in order to do this properly I need to ask Jeff to explain to me L254-259 in mri_prefers_3d/python/mri.py
Gs_local = np.zeros((len(kx_local), len(ts)), dtype=np.float64)
Gs_global = np.zeros((len(kx_global), len(ts)), dtype=np.float64)

nEVs_global = np.zeros_like(Gs_global)
nEVs_local = np.zeros_like(Gs_local)


def energy_norm_general(X1, X2, kx):
    u1 = X1['phi_z']
    w1 = -1.0j * kx * X1['phi']
    bx1 = X1['psi_z']
    bz1 = -1.0j * kx * X1['psi']
    u2 = X2['phi_z']
    w2 = -1.0j * kx * X2['phi']
    bx2 = X2['psi_z']
    bz2 = -1.0j * kx * X2['psi']
    field = (np.conj(u1)*u2 + np.conj(w1)*w2 + (1.0/MA2)*(np.conj(bx1)*bx2 + np.conj(bz1)*bz2)).evaluate().integrate()
    return field['g'][0]


z_basis = de.Chebyshev('z', Nz, interval=(-0.5 * Lz, 0.5 * Lz))
domain = de.Domain([z_basis], grid_dtype=np.complex128, comm=MPI.COMM_SELF)
probvars = ['phi', 'phi_z', 'psi', 'psi_z', 'uz', 'uzz']

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
problem.parameters['kx'] = kx0
problem.parameters['Re'] = Reynolds
problem.parameters['Rm'] = mReynolds
problem.substitutions['dx(A)'] = "1.0j*kx*A"
problem.substitutions['dt(A)'] = "gamma*A"
problem.substitutions['u'] = 'phi_z'
problem.substitutions['w'] = "-dx(phi)"
problem.substitutions['Bz'] = "-dx(psi)"
problem.substitutions['Bx'] = "psi_z"
problem.substitutions['Jz'] = "dx(dx(psi)) + dz(Bx)"
problem.substitutions['zeta'] = "dx(dx(phi)) + uz"
problem.substitutions['zetaz'] = "dx(dx(u)) + uzz"

problem.add_equation("dt(zeta) - 1/Re*(dx(dx(zeta)) + dz(zetaz)) + U*dx(zeta) - dx(phi)*dz(DU) - 1/MA2*(dx(Jz)) = 0")
problem.add_equation("uz - dz(u) = 0")
problem.add_equation("uzz - dz(uz) = 0")

problem.add_equation("dt(psi) - 1/Rm*Jz + U*dx(psi) - dx(phi) = 0")

problem.add_equation("phi_z - dz(phi) = 0")
problem.add_equation("psi_z - dz(psi) = 0")

problem.add_bc("left(w) = 0")
problem.add_bc("right(w) = 0")
problem.add_bc("left(dz(w)) = 0")  # no-slip, co-moving (with base flow U_0) boundaries
problem.add_bc("right(dz(w)) = 0")  # (does free-slip make more sense?)
problem.add_bc("right(Bz) = 0")  # perfectly conducting boundaries
problem.add_bc("left(Bz) = 0")

logger.info("done with BCs")

EP = Eigenproblem(problem, grow_func=lambda x: x.real, freq_func=lambda x: x.imag, reject=reject_spurious)
pencil = 0

for kx_ind, kx in enumerate(kx_local):
    logger.info("working on local calculation {}/{}".format(kx_ind+1, len(kx_local)))
    kxdir = kxdirs_local[kx_ind]
    try:
        Gs, eigenvalues, eigenvectors = LOP(EP, ts, {'kx': kx}, lambda x1, x2: energy_norm_general(x1, x2, kx), pencil=pencil, dense=dense, k=k)
    except ValueError:
        logger.warning("skipping kx = {} because no good eigenvalues were found".format(kx))
        Gs_local[kx_ind] = np.ones_like(Gs_local[kx_ind])*np.nan
        continue  # this means that no .h5 or .png files will be created for this ky, and the code will continue on to the next ky
    Gs_local[kx_ind] = Gs
    nEVs_local[kx_ind] = len(EP.evalues)
    with h5py.File(str(kxdir / "gains_mats.h5"), mode='w-') as file:
        setup_mats_grp = file.create_group('setup_mats')
        evalues_data = setup_mats_grp.create_dataset('eigenvalues', data=eigenvalues)

        scales_grp = file.create_group(
            'scales')  # in the spirit of "scales" in d2 h5 files, but this is a bad name for this use case
        ts_data = scales_grp.create_dataset('ts', data=ts)
        kx_data = scales_grp.create_dataset('kx', data=kx)
        Nz_data = scales_grp.create_dataset('Nz', data=Nz)
        Lz_data = scales_grp.create_dataset('Lz', data=Lz)
        zs_data = scales_grp.create_dataset('z', data=z)

        gains_grp = file.create_group('gains')
        Gs_data = gains_grp.create_dataset('Gs', data=Gs)
        # perturbations_data =
    if make_plots:
        plt.semilogy(ts, Gs)
        plt.semilogy(ts, Gs[0]*np.exp(ts*2.0*np.max(np.real(eigenvalues)))*Gs[-1]/np.exp(ts[-1]*2.0*np.max(np.real(eigenvalues))), '--')
        plt.xlabel(r'$T$')
        plt.ylabel(r'$G(T)$')
        plt.xlim((0, ts[-1]))
        plt.title(r'kx = {}'.format(kx))
        plt.savefig(str(kxdir / 'gain.png'))
        plt.close()

Gs_global[CW.rank::CW.size] = Gs_local
nEVs_global[CW.rank::CW.size] = nEVs_local
if CW.rank == 0:
    CW.Reduce(MPI.IN_PLACE, Gs_global, op=MPI.SUM, root=0)
    CW.Reduce(MPI.IN_PLACE, nEVs_global, op=MPI.SUM, root=0)
else:
    CW.Reduce(Gs_global, Gs_global, op=MPI.SUM, root=0)
    CW.Reduce(nEVs_global, nEVs_global, op=MPI.SUM, root=0)

if CW.rank == 0:
    with h5py.File(str(basedir / "gains.h5"), mode='w-') as file:
        gains_grp = file.create_group('gains')
        Gs_data = gains_grp.create_dataset('Gs', data=Gs_global)

        ks_grp = file.create_group('ks')
        kxs_data = ks_grp.create_dataset('kx', data=kx_global)

        scales_grp = file.create_group('scales')
        ts_data = scales_grp.create_dataset('ts', data=ts)
        Nz_data = scales_grp.create_dataset('Nz', data=Nz)
        Lz_data = scales_grp.create_dataset('Lz', data=Lz)
        zs_data = scales_grp.create_dataset('z', data=z)

        EVs_grp = file.create_group('EVs')
        nEVs_data = EVs_grp.create_dataset('nEVs', data=nEVs_global)
