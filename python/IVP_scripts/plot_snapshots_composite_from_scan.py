import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec as gridspec
import h5py
from configparser import ConfigParser
from pathlib import Path
from dedalus import public as de
import logging
logger = logging.getLogger(__name__)

run = "run_kx-lin-ky-lin_MA1_Re250_kmin0.1_Nz512-k600-saveLOPs"
LOP_t0_ind = 20
kxi = None
kyi = None

istep = 1500
field_names = ['u', 'v', 'w', 'Bx', 'By', 'Bz']

if kxi is None and kyi is None:
    get_max_ks = True
else:
    get_max_ks = False

# The .cfg file for this run
config_file = Path("../LOPs_3D/config_files/{}.cfg".format(run))
# data dir that was saved to
basedir = Path("../LOPs_3D") / config_file.stem
# Parse settings from .cfg file and save to variables
runconfig = ConfigParser()
runconfig.read(str(config_file))

# Parameters
params = runconfig['parameters']
Reynolds = params.getfloat('Reynolds')
MA = params.getfloat('MA')
MA2 = MA**2.0
Pm = params.getfloat('Pm')
mReynolds = Pm*Reynolds

t_start = params.getfloat('t_start')
t_stop = params.getfloat('t_stop')
t_num = params.getint('t_num')
t_endpoint = params.getboolean('t_endpoint')
times = np.linspace(t_start, t_stop, t_num, endpoint=t_endpoint)
save_LOPs_istep = params.getint('save_LOPs_istep')
LOP_ts = times[::save_LOPs_istep]

k = params.getint('k')
dense = params.getboolean('dense')
Nz = params.getint('Nz')
Lz = params.getfloat('Lz_factor')*np.pi

gains_file = "../LOPs_3D/{}/gains.h5".format(run)
with h5py.File(gains_file, mode='r') as file:
    kx_meshgrid = np.array(file['ks/kx'])
    ky_meshgrid = np.array(file['ks/ky'])
    kxs = kx_meshgrid[0]
    kys = ky_meshgrid[:, 0]
    Gs_global = np.array(file['gains/Gs'])
    Gs = Gs_global[:, :, ::save_LOPs_istep]
    # times = np.array(file['scales/ts'])

if get_max_ks:
    Gs_flat = Gs.reshape((-1, Gs.shape[-1]))  # flatten into a 2D array of shape (Gs.shape[0]*Gs.shape[1], Gs.shape[2])
    ind_flat = Gs_flat.argmax(axis=0)  # get the index where Gs_flat is maximized for each t
    kyind, kxind = np.unravel_index(ind_flat, Gs.shape[:2])  # unflatten/unravel those indices back into their kx, ky shape
    kxky_dirs = [basedir / 'kxi_{:04d}/kyi_{:04d}/'.format(kxind[i], kyind[i]) for i in range(len(Gs[0, 0]))]  # a list where each element points to the directory corresponding to the kx/ky that gives max G at each time
    kxky_dir = kxky_dirs[LOP_t0_ind]
    kxi = kxind[LOP_t0_ind]
    kyi = kyind[LOP_t0_ind]
    kx = kxs[kxi]
    ky = kys[kyi]
else:
    kxky_dir = basedir / 'kxi_{:04d}/kyi_{:04d}/'.format(kxi, kyi)
    kx = kxs[kxi]
    ky = kys[kyi]
out_dir = basedir / 'IVP_kxi{:04d}_kyi{:04d}_t0i{:03d}'.format(kxi, kyi, LOP_t0_ind)
plot_dir = out_dir / 'snapshots_composite'
if not plot_dir.exists():
    plot_dir.mkdir(parents=False)

with h5py.File(str(out_dir) + '/profiles_c/profiles_c_s1/profiles_c_s1_p0.h5', mode='r') as file:
    ts = np.array(file['scales/sim_time'])[::istep]
    field_data = dict([(field_names[i], np.array(file['tasks/{}'.format(field_names[i])])[::istep]) for i in range(len(field_names))])

z_basis = de.Chebyshev('z', Nz, interval=(-0.5 * Lz, 0.5 * Lz))
x_basis = de.Fourier('x', 4, interval=(0.0, 2.0 * np.pi / kx))
y_basis = de.Fourier('y', 4, interval=(0.0, 2.0 * np.pi / ky))
domain = de.Domain([x_basis, y_basis, z_basis], grid_dtype=np.float64)
# x = domain.grid(0)[:, 0, 0]
# y = domain.grid(1)[0, :, 0]
# z = domain.grid(2)[0, 0]
x1 = x_basis.grid(scale=1)
y1 = y_basis.grid(scale=1)
z1 = z_basis.grid(scale=1)
x16 = x_basis.grid(scale=16)
y16 = y_basis.grid(scale=16)
z16 = z_basis.grid(scale=16)

field = domain.new_field(name='field')
scale = 1.0

for tind, t in enumerate(ts):
    fig = plt.figure(figsize=(scale*12.0, 8.0))
    gs0 = gridspec.GridSpec(2, 3, figure=fig)
    gs0_subs = []
    # axes = []
    for fi, fname in enumerate(field_names):
        gs0_subs += [gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs0[fi])]
        # axes += [fig.add_subplot()]
        ax1 = fig.add_subplot(gs0_subs[-1][0, 0])
        ax2 = fig.add_subplot(gs0_subs[-1][0, 1])
        ax3 = fig.add_subplot(gs0_subs[-1][1, 0])
        ax4 = fig.add_subplot(gs0_subs[-1][1, 1])

        field.set_scales(1, keep_data=False)
        field['c'][1, 1, :] = field_data[fname][tind]
        field.set_scales(16, keep_data=True)
        field_out = field['g']

        pmesh1 = ax1.pcolormesh(x16, z16, field_out[:, 0, :].T)
        fig.colorbar(mappable=pmesh1, ax=ax1)
        ax1.set_title(fname+', $y=0$')
        ax1.tick_params(labelbottom=False)
        if fi == 0 or fi == 3:
            ax1.set_ylabel(r'$z$')

        pmesh2 = ax2.pcolormesh(y16, z16, field_out[0, :, :].T)
        fig.colorbar(mappable=pmesh2, ax=ax2)
        ax2.set_title(fname+', $x=0$')
        ax2.tick_params(labelbottom=False, labelleft=False)

        pmesh3 = ax3.pcolormesh(x16, z16, field_out[:, 15, :].T)
        fig.colorbar(mappable=pmesh3, ax=ax3)
        ax3.set_title(r'$y={:4.2f}$'.format(y16[15]))
        if fi == 0 or fi == 3:
            ax3.set_ylabel(r'$z$')
        else:
            ax1.tick_params(labelleft=False)
            ax3.tick_params(labelleft=False)
        if fi > 2:
            ax3.set_xlabel(r'$x$')

        pmesh4 = ax4.pcolormesh(y16, z16, field_out[15, :, :].T)
        fig.colorbar(mappable=pmesh4, ax=ax4)
        ax4.set_title(r'$x={:4.2f}$'.format(x16[15]))
        ax4.tick_params(labelleft=False)
        if fi > 2:
            ax4.set_xlabel(r'$y$')
        else:
            ax3.tick_params(labelbottom=False)
            ax4.tick_params(labelbottom=False)

    plt.suptitle(r'$t={:4.2f}$'.format(t))
    plt.tight_layout()
    plt.savefig(str(plot_dir) + '/snap_{:05d}'.format(tind))
    plt.close()
    if tind % 10 == 0:
        print(tind, "/", len(ts)-1)
