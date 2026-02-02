"""
Like plot_snapshots_composite_from_scan.py, except this makes quiver plots based on the u, w, bx, and bz components.
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec as gridspec
import h5py
from configparser import ConfigParser
from pathlib import Path
from dedalus import public as de
import logging
logger = logging.getLogger(__name__)

run = "run_kx-lin-ky-lin_MA1_Re250_kxmin0.1_Nz512-k600-saveLOPsAndEvolved"
LOP_t0_ind = 15
kxi = None
kyi = None
# fix_cbar = True

istep = 1
field_names = ['u', 'w', 'Bx', 'Bz']

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
t0 = LOP_ts[LOP_t0_ind]
# istop = np.argmin(np.abs(times-t0))


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
if get_max_ks:
    out_dir = basedir / 'IVP_kxi{:04d}_kyi{:04d}_t0i{:03d}'.format(kxi, kyi, LOP_t0_ind)
else:
    out_dir = basedir / 'IVP_specific_kxi{:04d}_kyi{:04d}_t0i{:03d}'.format(kxi, kyi, LOP_t0_ind)
# if fix_cbar:
#     plot_dir = out_dir / 'snapshots_composite_flat_fixed-cbar'
# else:
#     plot_dir = out_dir / 'snapshots_composite_flat'
plot_dir = out_dir / 'quiver_plots'
if not plot_dir.exists():
    plot_dir.mkdir(parents=False)

with h5py.File(str(out_dir) + '/profiles_c/profiles_c_s1/profiles_c_s1_p0.h5', mode='r') as file:
    full_ts = np.array(file['scales/sim_time'])
    istop = np.argmin(np.abs(full_ts - t0))
    ts = full_ts[:istop:istep]
    field_data = dict([(field_names[i], np.array(file['tasks/{}'.format(field_names[i])])[:istop:istep]) for i in range(len(field_names))])

z_basis = de.Chebyshev('z', Nz, interval=(-0.5 * Lz, 0.5 * Lz))
x_basis = de.Fourier('x', 32, interval=(0.0, 2.0 * np.pi / kx))
if ky > 0:
    y_basis = de.Fourier('y', 4, interval=(0.0, 2.0 * np.pi / ky))
    domain = de.Domain([x_basis, y_basis, z_basis], grid_dtype=np.float64)
    x1 = x_basis.grid(scale=1)
    y1 = y_basis.grid(scale=1)
    z1 = z_basis.grid(scale=1)
    # x16 = x_basis.grid(scale=16)
    # y16 = y_basis.grid(scale=16)
    # z16 = z_basis.grid(scale=16)
else:
    domain = de.Domain([x_basis, z_basis], grid_dtype=np.float64)
    x1 = x_basis.grid(scale=1)
    z1 = z_basis.grid(scale=1)
    # x16 = x_basis.grid(scale=16)
    # z16 = z_basis.grid(scale=16)
# x = domain.grid(0)[:, 0, 0]
# y = domain.grid(1)[0, :, 0]
# z = domain.grid(2)[0, 0]

field = domain.new_field(name='field')

field_data_g = {}
print(ts.shape)
for fi, fname in enumerate(field_names):
    print(fi)
    for tind, t in enumerate(ts):
        field.set_scales(1, keep_data=False)
        if ky > 0:
            field['c'][1, 1, :] = field_data[fname][tind]
        else:
            field['c'][1, :] = field_data[fname][tind]
        # field.set_scales(16, keep_data=True)
        field.set_scales(1, keep_data=True)
        # field_out = field['g']
        if tind == 0:
            if ky > 0:
                field_data_g[fname] = np.zeros((len(ts), *field['g'][:, 0, :].shape), dtype=np.float64)
            else:
                field_data_g[fname] = np.zeros((len(ts), *field['g'].shape), dtype=np.float64)
        if ky > 0:
            field_data_g[fname][tind] = field['g'][:, 0, :]
        else:
            field_data_g[fname][tind] = field['g']


scale = 1.0

for tind, t in enumerate(ts):
    fig = plt.figure(figsize=(scale*8.0, 8.0))
    # gs0 = gridspec.GridSpec(2, 1, figure=fig)
    # gs0_subs = []
    # # axes = []
    # for fi, fname in enumerate(field_names):
    #     gs0_subs += [gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[fi])]
    #     ax1 = fig.add_subplot(gs0_subs[-1][0, 0])
    #     # if fix_cbar:
    #     #     field_max = np.max(np.abs(field_data_g[fname][:, :-1, :-1]))
    #     #     pmesh1 = ax1.pcolormesh(x16, z16, field_data_g[fname][tind][:-1, :-1].T, shading='flat', vmin=-field_max, vmax=field_max)
    #     # else:
    #     #     pmesh1 = ax1.pcolormesh(x16, z16, field_data_g[fname][tind][:-1, :-1].T, shading='flat')
    #     # fig.colorbar(mappable=pmesh1, ax=ax1)
    #
    #     ax1.set_title(fname+', $y=0$')
    #     if fi == 0 or fi == 3:
    #         ax1.set_ylabel(r'$z$')
    #     ax1.set_xlabel(r'$x$')
    # plt.suptitle(r'$t={:4.2f}$'.format(t))
    plt.subplot(2, 1, 1)
    # ax = plt.gca()
    zskip = 16
    zstart = 150
    zstop = -150
    plt.quiver(x1, z1[zstart:zstop:zskip], field_data_g['u'][tind][:, zstart:zstop:zskip].T, field_data_g['w'][tind][:, zstart:zstop:zskip].T)#, scale=10.0)
    plt.subplot(2, 1, 2)
    plt.quiver(x1, z1[zstart:zstop:zskip], field_data_g['Bx'][tind][:, zstart:zstop:zskip].T, field_data_g['Bz'][tind][:, zstart:zstop:zskip].T)#, scale=10.0)

    plt.tight_layout()
    plt.savefig(str(plot_dir) + '/snap_{:05d}'.format(tind))
    plt.close()
    if tind % 10 == 0:
        print(tind, "/", len(ts)-1)
