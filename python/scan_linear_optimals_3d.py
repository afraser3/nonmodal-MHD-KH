"""
Does all the stuff in LOP_save_mats but as a kx, ky scan.

To run in parallel with 2 processes, e.g., do:
    $ mpiexec -n 2 python3 scan_linear_optimals_3d.py LOPs_3D/config_files/run1.cfg
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
# data dir to save to: if "LOPs_3D/config_files/run1.cfg" is input, then "LOPs_3D/run1/" becomes the output directory
basedir = Path("LOPs_3D") / config_file.stem
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
kx_log = params.getboolean('kx_log')
kx_start = params.getfloat('kx_start')
kx_stop = params.getfloat('kx_stop')
kx_num = params.getint('kx_num')
kx_endpoint = params.getboolean('kx_endpoint')
add_kx_zero = params.getboolean('add_kx_zero')  # this is ugly but a lazy way to add kx=0 to log scans!
ky_start = params.getfloat('ky_start')
ky_stop = params.getfloat('ky_stop')
ky_num = params.getint('ky_num')
ky_endpoint = params.getboolean('ky_endpoint')

t_start = params.getfloat('t_start')
t_stop = params.getfloat('t_stop')
t_num = params.getint('t_num')
t_endpoint = params.getboolean('t_endpoint')
ts = np.linspace(t_start, t_stop, t_num, endpoint=t_endpoint)

k = params.getint('k')
dense = params.getboolean('dense')
Nz = params.getint('Nz')
Lz = params.getfloat('Lz_factor')*np.pi
Legendre = params.getboolean('Legendre')
reject_spurious = params.getboolean('reject_spurious')
make_plots = params.getboolean('make_plots')
plot_modes = params.getboolean('plot_modes')
save_evecs = params.getboolean('save_evecs')
save_LOPs = params.getboolean('save_LOPs')
save_LOPs_istep = params.getint('save_LOPs_istep')
save_betas = params.getboolean('save_betas')  # the amplitude of the LOPs expressed in the eigenmode basis

# kx = 0.2
# ky = 1e-8  # worked for ky=0.1, but 0.2 gives no good eigenvalues, strangely. At Nz=1024 it gave 1 good eigenvalue.
# MA = 10.0  # Alfven Mach number. For this value, the system is unstable.

ky_global = np.linspace(ky_start, ky_stop, ky_num, endpoint=ky_endpoint)  # the full set of kys to scan over
if kx_log:
    kx_global = np.geomspace(kx_start, kx_stop, kx_num, endpoint=kx_endpoint)  # the full set of kxs to scan over
else:
    kx_global = np.linspace(kx_start, kx_stop, kx_num, endpoint=kx_endpoint)  # the full set of kxs to scan over
if add_kx_zero and kx_global[0] != 0:
    kx_global = np.insert(kx_global, 0, 0.0)
kxdirs_global = [basedir / 'kxi_{:04d}/'.format(kx_i) for kx_i, kx in enumerate(kx_global)]
if CW.rank == 0:
    for kxdir in kxdirs_global:
        if not kxdir.exists():
            kxdir.mkdir(parents=True)

kx0 = kx_global[0]  # just a placeholder MA for setting Dedalus EVP up
ky0 = ky_global[0]
if CW.rank == 0:
    logger.info('Scanning over these kxs: {}'.format(kx_global))
    logger.info('And scanning over these kys: {}'.format(ky_global))

# kx_local = kx_global[CW.rank::CW.size]  # the kxs that this processor will scan over
# logger.debug('This MPI process will scan over these kxs: {}'.format(kx_local))

# in order to do this properly I need to ask Jeff to explain to me L254-259 in mri_prefers_3d/python/mri.py
# Gs_local = np.zeros_like((len(kx_local), len(ky_global), len(ts)), dtype=np.float64)
# Gs_global = np.zeros_like((len(kx_global), len(ky_global), len(ts)), dtype=np.float64)


def energy_norm_general(X1, X2, kx, ky):
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

    field = (np.conj(u1)*u2 + np.conj(v1)*v2 + np.conj(w1)*w2 + (1.0/MA2)*(np.conj(bx1)*bx2 + np.conj(by1)*by2 + np.conj(bz1)*bz2)).evaluate().integrate()
    return field['g'][0]


if Legendre:
    z_basis = de.Legendre('z', Nz, interval=(-0.5 * Lz, 0.5 * Lz))
else:
    z_basis = de.Chebyshev('z', Nz, interval=(-0.5 * Lz, 0.5 * Lz))
domain = de.Domain([z_basis], grid_dtype=np.complex128, comm=MPI.COMM_SELF)
probvars = ['u', 'v', 'w', 'p', 'Ax', 'Ay', 'Az', 'Phi', 'omega_x', 'omega_y', 'Bx', 'By']

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
problem.parameters['ky'] = ky0
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

parameter_grids = np.meshgrid(kx_global, ky_global)
dims = parameter_grids[0].shape
index = np.arange(np.prod(dims))
load_indices = np.array_split(index, CW.size)
my_indices = load_indices[CW.rank]

local_Gs = np.zeros((my_indices.size, len(ts)), dtype=np.float64)
local_kx = np.zeros(my_indices.size, dtype=np.float64)
local_ky = np.zeros(my_indices.size, dtype=np.float64)
local_nEVs = np.zeros(my_indices.size, dtype=np.float64)  # should really be int but if it ain't broke don't fix it
local_cond = np.zeros_like(local_nEVs)
for n, index in enumerate(my_indices):
    logger.info("working on local calculation {}/{}".format(n+1, len(my_indices)))
    unraveled_index = np.unravel_index(index, dims)
    values = [parameter_grids[i][unraveled_index] for i, v in enumerate(parameter_grids)]
    kx = values[0]
    local_kx[n] = kx
    kxi = unraveled_index[1]
    ky = values[1]
    local_ky[n] = ky
    kyi = unraveled_index[0]
    local_dir = basedir / 'kxi_{:04d}'.format(kxi) / 'kyi_{:04d}'.format(kyi)
    if kx > 0. or ky > 0.:
        try:
            # Gs, eigenvalues, eigenvectors, F, Q, R, pre_right = LOP(EP, ts, {'kx': kx, 'ky': ky}, lambda x1, x2: energy_norm_general(x1, x2, kx, ky), pencil=pencil, dense=dense, k=k, return_setup_mats=True)
            # Gs, perturbations, eigenvalues, eigenvectors = LOP(EP, ts, {'kx': kx, 'ky': ky}, lambda x1, x2: energy_norm_general(x1, x2, kx, ky), pencil=pencil, dense=dense, k=k, return_perturbations=True)
            if save_LOPs:
                if save_betas:
                    Gs, perturbations, eigenvalues, eigenvectors = LOP(EP, ts, {'kx': kx, 'ky': ky}, lambda x1, x2: energy_norm_general(x1, x2, kx, ky), pencil=pencil, dense=dense, k=k, return_perturbations=True, left=True, normalize_left=False)
                else:
                    Gs, perturbations, eigenvalues, eigenvectors = LOP(EP, ts, {'kx': kx, 'ky': ky}, lambda x1, x2: energy_norm_general(x1, x2, kx, ky), pencil=pencil, dense=dense, k=k, return_perturbations=True)
            else:
                Gs, eigenvalues, eigenvectors = LOP(EP, ts, {'kx': kx, 'ky': ky}, lambda x1, x2: energy_norm_general(x1, x2, kx, ky), pencil=pencil, dense=dense, k=k)
        except ValueError:
            logger.warning("skipping (kx, ky) = ({}, {}) because no good eigenvalues were found".format(kx, ky))
            local_Gs[n] = np.ones_like(local_Gs[n])*np.nan
            if not local_dir.exists():
                local_dir.mkdir(parents=True)
            with h5py.File(str(local_dir / "gains_mats.h5"), mode='w-') as file:
                setup_mats_grp = file.create_group('setup_mats')
                evalues_data = setup_mats_grp.create_dataset('eigenvalues', data=EP.evalues)
                evalues_high_data = setup_mats_grp.create_dataset('eigenvalues_high', data=EP.evalues_high)
                evalues_low_data = setup_mats_grp.create_dataset('eigenvalues_low', data=EP.evalues_low)
                if save_evecs:
                    evecs_data = setup_mats_grp.create_dataset('eigenvectors', data=((EP.solver.eigenvectors.T)[EP.evalues_index]).T)
                scales_grp = file.create_group('scales')
                ts_data = scales_grp.create_dataset('ts', data=ts)
                kx_data = scales_grp.create_dataset('kx', data=kx)
                ky_data = scales_grp.create_dataset('ky', data=ky)
                Nz_data = scales_grp.create_dataset('Nz', data=Nz)
                Lz_data = scales_grp.create_dataset('Lz', data=Lz)
                zs_data = scales_grp.create_dataset('z', data=z)
            continue  # this means that no .h5 or .png files will be created for this ky, and the code will continue on to the next ky
        except np.linalg.LinAlgError:
            # If the Cholesky decomposition doesn't work, say so, and dump some data to help diagnose whatever is going on
            logger.warning("skipping (kx, ky) = ({}, {}) because LinAlgError was returned -- look into this!".format(kx, ky))
            if not local_dir.exists():
                local_dir.mkdir(parents=True)
            with h5py.File(str(local_dir / "gains_mats.h5"), mode='w-') as file:
                setup_mats_grp = file.create_group('setup_mats')
                evalues_data = setup_mats_grp.create_dataset('eigenvalues', data=EP.evalues)
                evalues_high_data = setup_mats_grp.create_dataset('eigenvalues_high', data=EP.evalues_high)
                evalues_low_data = setup_mats_grp.create_dataset('eigenvalues_low', data=EP.evalues_low)
                if save_evecs:
                    evecs_data = setup_mats_grp.create_dataset('eigenvectors', data=((EP.solver.eigenvectors.T)[EP.evalues_index]).T)
                scales_grp = file.create_group('scales')
                ts_data = scales_grp.create_dataset('ts', data=ts)
                kx_data = scales_grp.create_dataset('kx', data=kx)
                ky_data = scales_grp.create_dataset('ky', data=ky)
                Nz_data = scales_grp.create_dataset('Nz', data=Nz)
                Lz_data = scales_grp.create_dataset('Lz', data=Lz)
                zs_data = scales_grp.create_dataset('z', data=z)
            continue
    else:
        local_Gs[n] = np.ones_like(local_Gs[n])*np.nan
        continue
    # Gs_local[kx_ind, ky_ind] = Gs
    local_Gs[n] = Gs
    local_nEVs[n] = len(EP.evalues)
    local_cond[n] = np.linalg.cond(EP.solver.pencils[pencil].L_exp.toarray())
    if not local_dir.exists():
        local_dir.mkdir(parents=True)
    with h5py.File(str(local_dir / "gains_mats.h5"), mode='w-') as file:
        setup_mats_grp = file.create_group('setup_mats')
        evalues_data = setup_mats_grp.create_dataset('eigenvalues', data=eigenvalues)
        evalues_high_data = setup_mats_grp.create_dataset('eigenvalues_high', data=EP.evalues_high)
        evals_low = EP.evalues_low
        if dense:  # this low_order stuff doesn't matter if using dense solve because lefts and rights are in the same order
            low_order = np.array(range(len(evals_low)))
        else:
            low_order = np.argsort(np.imag(evals_low))  # reorder low evals so they match the ordering of the betas, below
        evalues_low_data = setup_mats_grp.create_dataset('eigenvalues_low', data=evals_low[low_order])
        if save_evecs:
            evecs_data = setup_mats_grp.create_dataset('eigenvectors', data=eigenvectors)
        if save_betas:
            # left_order = np.argsort(-np.imag(EP.solver.left_eigenvalues))  # reorder left evals so they match the ordering of the betas, below
            if dense:
                left_order = low_order
            else:
                left_order = np.zeros_like(low_order)
                left_evals = EP.solver.left_eigenvalues
                for i in range(len(left_order)):
                    left_order[i] = np.argmin(np.abs(np.conj(left_evals) - evals_low[low_order[i]]))
                left_evalues_data = setup_mats_grp.create_dataset('left_eigenvalues', data=left_evals[left_order])
            evecs_low_data = setup_mats_grp.create_dataset('eigenvectors_low', data=EP.solver.eigenvectors[:, low_order])  # just for troubleshooting for now
        # R_data = setup_mats_grp.create_dataset('R', data=R)
        # Q_data = setup_mats_grp.create_dataset('Q', data=Q)
        # F_data = setup_mats_grp.create_dataset('F', data=F)
        # pre_right_data = setup_mats_grp.create_dataset('pre_right', data=pre_right.toarray())

        scales_grp = file.create_group('scales')  # in the spirit of "scales" in d2 h5 files, but this is a bad name for this use case
        ts_data = scales_grp.create_dataset('ts', data=ts)
        kx_data = scales_grp.create_dataset('kx', data=kx)
        ky_data = scales_grp.create_dataset('ky', data=ky)
        Nz_data = scales_grp.create_dataset('Nz', data=Nz)
        Lz_data = scales_grp.create_dataset('Lz', data=Lz)
        zs_data = scales_grp.create_dataset('z', data=z)

        gains_grp = file.create_group('gains')
        Gs_data = gains_grp.create_dataset('Gs', data=Gs)
        if save_LOPs:
            LOPs_data = gains_grp.create_dataset('LOPs', data=perturbations[:, ::save_LOPs_istep])
            if save_betas:
                modified_lefts_data = gains_grp.create_dataset('modified_left_eigenvectors', data=EP.solver.modified_left_eigenvectors)  # just for troubleshooting for now
                modified_lefts = EP.solver.modified_left_eigenvectors
                # Ideally would like to normalize the modified left eigenvectors such that they are orthoNORMAL, not
                # just orthogonal, to the right eigenvectors. Doing that requires first making sure that the left and
                # right eigenvectors are arranged in the same order. This assumes that the left and right eigenvalues
                # are indeed the same just not arranged in the same order -- might be worth raising a warning if they
                # aren't actually the same
                evecs_low_reorder = EP.solver.eigenvectors[:, low_order]  # reordered right eigenvectors
                modified_left_evecs_reorder = modified_lefts[:, left_order]  # reordered modified left evecs
                # next, normalize the left evecs
                modified_left_evecs_reorder_norm = modified_left_evecs_reorder / np.conjugate(np.diag(modified_left_evecs_reorder.T.conj() @ evecs_low_reorder))
                # betas = np.matmul(np.conjugate(EP.solver.modified_left_eigenvectors.T), perturbations[:, ::save_LOPs_istep])
                betas = modified_left_evecs_reorder_norm.T.conj() @ perturbations[:, ::save_LOPs_istep]
                betas_data = gains_grp.create_dataset('betas', data=betas)
    if make_plots:
        plt.semilogy(ts, Gs)
        plt.semilogy(ts, Gs[0]*np.exp(ts*2.0*np.max(np.real(eigenvalues)))*Gs[-1]/np.exp(ts[-1]*2.0*np.max(np.real(eigenvalues))), '--')
        plt.xlabel(r'$T$')
        plt.ylabel(r'$G(T)$')
        plt.xlim((0, ts[-1]))
        plt.title(r'kx = {}, ky = {}'.format(kx, ky))
        plt.savefig(str(local_dir / 'gain.png'))
        plt.close()
    if plot_modes:
        emode_dir = local_dir / 'emode_plots'
        if not emode_dir.exists():
            emode_dir.mkdir(parents=True)
        for ev_ind, ev in enumerate(eigenvalues):
            fig = EP.plot_mode(ev_ind)
            fig.suptitle('EV: {}'.format(ev))
            fig.savefig(str(emode_dir / 'emode_{}.pdf'.format(ev_ind)))
            plt.close()

gathered_Gs = np.zeros((*dims, len(ts)), dtype=np.float64)
rec_counts = np.array([s.size for s in load_indices])*len(ts)
displacements = np.cumsum(rec_counts) - rec_counts
CW.Gatherv(local_Gs, [gathered_Gs, rec_counts, displacements, MPI.F_DOUBLE])

rec_counts_kx = np.array([s.size for s in load_indices])
displacements_kx = np.cumsum(rec_counts_kx) - rec_counts_kx
gathered_kx = np.zeros(dims, dtype=np.float64)
gathered_ky = np.zeros(dims, dtype=np.float64)
gathered_nEVs = np.zeros(dims, dtype=np.float64)
gathered_cond = np.zeros(dims, dtype=np.float64)
CW.Gatherv(local_kx, [gathered_kx, rec_counts_kx, displacements_kx, MPI.F_DOUBLE])
CW.Gatherv(local_ky, [gathered_ky, rec_counts_kx, displacements_kx, MPI.F_DOUBLE])
CW.Gatherv(local_nEVs, [gathered_nEVs, rec_counts_kx, displacements_kx, MPI.F_DOUBLE])
CW.Gatherv(local_cond, [gathered_cond, rec_counts_kx, displacements_kx, MPI.F_DOUBLE])
# print(CW.rank, rec_counts, displacements)
# if CW.rank == 0:
    # print(np.shape(data))
    # print(np.shape(local_Gs))
# print(CW.rank, local_Gs)
if CW.rank == 0:
    with h5py.File(str(basedir / "gains.h5"), mode='w-') as file:
        gains_grp = file.create_group('gains')
        Gs_data = gains_grp.create_dataset('Gs', data=gathered_Gs)

        ks_grp = file.create_group('ks')
        kxs_data = ks_grp.create_dataset('kx', data=gathered_kx)
        kys_data = ks_grp.create_dataset('ky', data=gathered_ky)

        scales_grp = file.create_group('scales')
        ts_data = scales_grp.create_dataset('ts', data=ts)
        Nz_data = scales_grp.create_dataset('Nz', data=Nz)
        Lz_data = scales_grp.create_dataset('Lz', data=Lz)
        zs_data = scales_grp.create_dataset('z', data=z)

        EVs_grp = file.create_group('EVs')
        nEVs_data = EVs_grp.create_dataset('nEVs', data=gathered_nEVs)
        cond_data = EVs_grp.create_dataset('cond', data=gathered_cond)

# fix data saving scheme (ditch the individual .h5's, and include number of good eigenmodes and condition number in saved data)
# save nans on bad evalues
# plot number of eigenvalues
# plot condition number
# scan from k_x,y = 0.1 to 0.2 to see if we can resolve a ridge
# for 0.25, 0.26 (where we get 100 good eigenvalues) are the good eigenmodes really sharp? do they get sharper as ky -> kx?
# in plot_mode, there's a phase-setter that sets a field to real
# put [specifyable_index] at the end of line 293 in eigenproblem.py
# check if the block structure moves around if MA changes very slightly
