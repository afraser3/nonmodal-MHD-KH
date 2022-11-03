import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec as gridspec
import h5py
from pathlib import Path
from dedalus import public as de
import logging
logger = logging.getLogger(__name__)

test_ind = 31

Nz = 512
Reynolds = 250.0
Pm = 1.0
mReynolds = Pm * Reynolds
Lz = 10.0*np.pi
kx = 0.025
ky = 0.6
MA = 0.25
MA2 = MA**2.0
k = 300
dense = False

field_names = ['u', 'v', 'w', 'Bx', 'By', 'Bz']
istep = 1500

if dense:
    h5dir = Path("LOPs_3D/individual/Re{}/MA{}/Pm{}/kx{}/ky{}".format(Reynolds, MA, Pm, kx, ky))
    h5_fname = str(h5dir / 'Nz{}_dense.h5'.format(Nz))
    out_dir = h5dir / 'Nz{}_dense_ti{}'.format(Nz, test_ind)
    plot_dir = out_dir / 'snapshots_composite'
else:
    h5dir = Path("LOPs_3D/individual/Re{}/MA{}/Pm{}/kx{}/ky{}".format(Reynolds, MA, Pm, kx, ky))
    h5_fname = str(h5dir / 'Nz{}_k{}.h5'.format(Nz, k))
    out_dir = h5dir / 'Nz{}_k{}_ti{}'.format(Nz, k, test_ind)
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
