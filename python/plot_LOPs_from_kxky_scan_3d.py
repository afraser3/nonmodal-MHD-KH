# Reads the files made by scan_linear_optimals_3d.py and plots the perturbations
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import h5py
from dedalus import public as de
from pathlib import Path
from configparser import ConfigParser
import glob
import re


def natural_sort(iterable, reverse=False):
    """
    Sort alphanumeric strings naturally, i.e. with "1" before "10".
    Based on http://stackoverflow.com/a/4836734.

    Actually, just copy+pasted from dedalus
    """

    convert = lambda sub: int(sub) if sub.isdigit() else sub.lower()
    key = lambda item: [convert(sub) for sub in re.split('([0-9]+)', str(item))]

    return sorted(iterable, key=key, reverse=reverse)


# run = 'run_kx-log-scan-saveLOPs-saveEvecs-saveBetas'
# run = "run_kx-lin-ky-lin_MA1_Re250_kxmin0.1_Nz512-k600-saveLOPsAndEvolved"
# run = "run_kx-lin-ky-lin_MA1.1_Re250_kxmin0.1_T200_Nz512-k600-saveLOPsAndEvolved"
run = "run_kx-lin-ky-lin_MA0.5_Re250_kxmin0.1_T100_Nz512-k600-saveLOPsAndEvolved"

sparse = False
basedir = Path("LOPs_3D") / run
config_file = Path("LOPs_3D/config_files/"+run+'.cfg')
# gains_files = glob.glob(str(basedir)+'/kxi*/kyi*/gains_mats.h5')
# kydirs = glob.glob(str(basedir)+'/kxi*/kyi*/')

# Parse settings from .cfg file and save to variables
runconfig = ConfigParser()
runconfig.read(str(config_file))

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
ky_start = params.getfloat('ky_start')
ky_stop = params.getfloat('ky_stop')
ky_num = params.getint('ky_num')
ky_endpoint = params.getboolean('ky_endpoint')
kys = np.linspace(ky_start, ky_stop, ky_num, endpoint=ky_endpoint)  # the full set of kys to scan over
# TODO: need to figure out what to do if kxs[0] == 0, since the rest of this script assumes finite kx
# I guess the same goes for ky
if kx_log:
    kxs = np.geomspace(kx_start, kx_stop, kx_num, endpoint=kx_endpoint)  # the full set of kxs to scan over
else:
    kxs = np.linspace(kx_start, kx_stop, kx_num, endpoint=kx_endpoint)  # the full set of kxs to scan over
# kxdirs = [basedir / 'kxi_{:04d}/'.format(kx_i) for kx_i, kx in enumerate(kxs)]
# kydirs = [Path(kydir) for kydir in natural_sort(glob.glob('LOPs_3D/{}/kxi_*/kyi_0000'.format(run)))]

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
make_plots = params.getboolean('make_plots')
plot_modes = params.getboolean('plot_modes')
save_betas = params.getboolean('save_betas')
save_LOPs = params.getboolean('save_LOPs')
save_LOPs_istep = params.getint('save_LOPs_istep')
LOP_ts = ts[::save_LOPs_istep]

global_gains_file = "LOPs_3D/{}/gains.h5".format(run)
with h5py.File(global_gains_file, mode='r') as file:
    # kxs = np.array(file['ks/kx'])
    # kys = np.array(file['ks/ky'])
    Gs = np.array(file['gains/Gs'])[:, :, ::save_LOPs_istep]
    # times = np.array(file['scales/ts'])
Gs_flat = Gs.reshape((-1, Gs.shape[-1]))  # flatten into a 2D array of shape (Gs.shape[0]*Gs.shape[1], Gs.shape[2])
ind_flat = Gs_flat.argmax(axis=0)  # get the index where Gs_flat is maximized for each t
kyind, kxind = np.unravel_index(ind_flat, Gs.shape[:2])  # unflatten/unravel those indices back into their kx, ky shape
kxky_dirs = [basedir / 'kxi_{:04d}/kyi_{:04d}/'.format(kxind[i], kyind[i]) for i in range(len(Gs[0, 0]))]  # a list where each element points to the directory corresponding to the kx/ky that gives max G at each time
# tind = 2
# target_time = LOP_ts[tind]
if params.getboolean('dense'):
    solver = 'dense'
else:
    solver = 'k{}'.format(params['k'])
if params.getboolean('Legendre'):
    basis = 'Legendre'
else:
    basis = 'Chebyshev'
# plotname = 'LOPs-vs-t0_MA{}_Nz{}_{}_{}.pdf'.format(params['MA'], params['Nz'], solver, basis)
plotname = 'LOPs-vs-t0__{}.pdf'.format(run)
LOP_dir = basedir/'LOPs'
if not LOP_dir.exists():
    LOP_dir.mkdir(parents=False)

with PdfPages(str(basedir) + '/' + plotname) as pdf:
    # for kydir in kydirs[-2::-1]:  # for this one case, skip this one ky because it's not behaving well
    # for kydir in kydirs[::-1]:
    for ti, kxkydir in enumerate(kxky_dirs):
        if LOP_ts[ti] == 0:  # no point in plotting anything for t=0
            continue
        print('starting: ', kxkydir)
        gains_file = str(kxkydir / 'gains_mats.h5')
        with h5py.File(gains_file, mode='r') as file:
            kx = np.array(file['scales/kx'])
            ky = np.array(file['scales/ky'])
            try:
                LOPs = np.array(file['gains/LOPs'])
            except KeyError:
                continue
            if sparse and save_betas:
                left_evals = np.array(file['setup_mats/left_eigenvalues'])
            evals = np.array(file['setup_mats/eigenvalues'])
            evals_low = np.array(file['setup_mats/eigenvalues_low'])
            if save_betas:
                betas = np.array(file['gains/betas'])
            # loop over the t axis of the LOPs array and shift the complex phase of the LOPs
            # (this can certainly be vectorized)
            # for ti in range(np.shape(LOPs)[-1]):
                # LOP_phase = LOPs[0, ti]/np.abs(LOPs[0, ti])
                # LOPs[:, ti] = LOPs[:, ti]/LOP_phase

        LOP = LOPs[:, ti]  # each kx/ky includes the LOPs at several times, but I just want to grab the LOP corresponding the time at which this kx/ky is where G maximizes
        if save_betas:
            betas = betas[:, ti]

        if Legendre:
            z_basis = de.Legendre('z', Nz, interval=(-0.5 * Lz, 0.5 * Lz))
        else:
            z_basis = de.Chebyshev('z', Nz, interval=(-0.5 * Lz, 0.5 * Lz))
        x_basis = de.Fourier('x', 32, interval=(0.0, 2.0 * np.pi / kx))
        if ky > 0:
            y_basis = de.Fourier('x', 32, interval=(0.0, 2.0 * np.pi / ky))
        else:
            y_basis = de.Fourier('x', 4, interval=(0.0, 2.0 * np.pi / kx))
        domain = de.Domain([x_basis, y_basis, z_basis], grid_dtype=np.float64)
        x = domain.grid(0)[:, 0, 0]
        y = domain.grid(1)[0, :, 0]
        z = domain.grid(2)[0, 0, :]
        probvars = ['u', 'v', 'w', 'p', 'Ax', 'Ay', 'Az', 'Phi', 'omega_x', 'omega_y', 'Bx', 'By']
        LOP_fields = [domain.new_field(name=probvar) for probvar in probvars]
        for ind, field in enumerate(LOP_fields):
            # select the part of the array that corresponds to each field
            slice_start = int(ind*len(LOP)/len(probvars))
            slice_stop = int((ind+1)*len(LOP)/len(probvars))
            if ky > 0:
                field['c'][1, 1, :] = LOP[slice_start:slice_stop]
            else:
                field['c'][1, 0, :] = LOP[slice_start:slice_stop]
        Bz = domain.new_field(name='Bz')
        Bz['c'] = 1j*kx*LOP_fields[5]['c'] - 1j*ky*LOP_fields[4]['c']  # Bz = dx(Ay) - dy(Ax)

        yind = 0  # where in y to slice
        scale = 0.8
        fig = plt.figure(figsize=(scale*16, scale*10))

        plt.subplot(2, 4, 1)
        plt.pcolormesh(x, z, LOP_fields[0]['g'][:, yind, :].T, rasterized=True)
        # plt.title(r'$u$, $t={}$'.format(LOP_ts[tind]))
        plt.title(r'$u$')
        plt.ylabel(r'$z$')
        plt.colorbar()

        plt.subplot(2, 4, 2)
        plt.pcolormesh(x, z, LOP_fields[1]['g'][:, yind, :].T, rasterized=True)
        plt.title(r'$v$')
        # plt.ylabel(r'$z$')
        plt.colorbar()

        plt.subplot(2, 4, 3)
        plt.pcolormesh(x, z, LOP_fields[2]['g'][:, yind, :].T, rasterized=True)
        plt.title(r'$w$')
        # plt.ylabel(r'$z$')
        plt.colorbar()

        plt.subplot(2, 4, 5)
        plt.pcolormesh(x, z, LOP_fields[10]['g'][:, yind, :].T, rasterized=True)
        plt.title(r'$B_x$')
        plt.ylabel(r'$z$')
        plt.xlabel(r'$x$')
        plt.colorbar()

        plt.subplot(2, 4, 6)
        plt.pcolormesh(x, z, LOP_fields[11]['g'][:, yind, :].T, rasterized=True)
        plt.title(r'$B_y$')
        # plt.ylabel(r'$z$')
        plt.xlabel(r'$x$')
        plt.colorbar()

        plt.subplot(2, 4, 7)
        plt.pcolormesh(x, z, Bz['g'][:, yind, :].T, rasterized=True)
        plt.title(r'$B_z$')
        # plt.ylabel(r'$z$')
        plt.xlabel(r'$x$')
        plt.colorbar()

        if save_betas:
            plt.subplot(2, 4, 4)
            # plt.plot(np.real(evals), np.imag(evals), 'x')
            beta_order = np.argsort(np.abs(betas))
            if sparse:
                plt.scatter(np.real(left_evals[beta_order]), -np.imag(left_evals[beta_order]), marker='.', c=np.abs(betas[beta_order]))
            else:
                plt.scatter(np.real(evals_low[beta_order]), -np.imag(evals_low[beta_order]), marker='.', c=np.abs(betas[beta_order]))
            plt.colorbar()

        plt.suptitle(r'$k_x={:1.2e}$, $k_y={:1.2e}$, $t_0={}$'.format(kx, ky, LOP_ts[ti]))

        plt.tight_layout()
        # plt.show()
        # plt.savefig(str(kydir / 'LOP_tind{}.png'.format(tind)))
        # plt.savefig(str(kydir / 'LOP_tind{}.pdf'.format(tind)))
        pdf.savefig()
        plt.savefig(str(LOP_dir) + '/LOP_ti{}__{}.pdf'.format(ti, run))
        plt.savefig(str(LOP_dir) + '/LOP_ti{}__{}.png'.format(ti, run))
        plt.close()
