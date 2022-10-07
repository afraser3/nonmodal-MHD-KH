import numpy as np
from matplotlib import pyplot as plt
import h5py
from pathlib import Path
from dedalus import public as de
import logging
logger = logging.getLogger(__name__)

test_ind = 50

Nz = 256
Reynolds = 50.0
Pm = 1.0
mReynolds = Pm * Reynolds
Lz = 10.0*np.pi
kx = 0.2  # 0.025
ky = 0.025  # 0.35
MA = 1.0
MA2 = MA**2.0
k = 300
dense = False

field1_name = 'w'
field2_name = 'Bz'
istep = 2

if dense:
    h5dir = Path("LOPs_3D/individual/Re{}/MA{}/Pm{}/kx{}/ky{}".format(Reynolds, MA, Pm, kx, ky))
    h5_fname = str(h5dir / 'Nz{}_dense.h5'.format(Nz))
    out_dir = h5dir / 'Nz{}_dense_ti{}'.format(Nz, test_ind)
    plot_dir = out_dir / 'snapshots'
else:
    h5dir = Path("LOPs_3D/individual/Re{}/MA{}/Pm{}/kx{}/ky{}".format(Reynolds, MA, Pm, kx, ky))
    h5_fname = str(h5dir / 'Nz{}_k{}.h5'.format(Nz, k))
    out_dir = h5dir / 'Nz{}_k{}_ti{}'.format(Nz, k, test_ind)
    plot_dir = out_dir / 'snapshots'
if not plot_dir.exists():
    plot_dir.mkdir(parents=False)

with h5py.File(str(out_dir) + '/profiles_c/profiles_c_s1/profiles_c_s1_p0.h5', mode='r') as file:
    ts = np.array(file['scales/sim_time'])[::istep]
    field1_data = np.array(file['tasks/{}'.format(field1_name)])[::istep]
    field2_data = np.array(file['tasks/{}'.format(field1_name)])[::istep]

z_basis = de.Chebyshev('z', Nz, interval=(-0.5 * Lz, 0.5 * Lz))
x_basis = de.Fourier('x', 32, interval=(0.0, 2.0 * np.pi / kx))
y_basis = de.Fourier('y', 32, interval=(0.0, 2.0 * np.pi / ky))
domain = de.Domain([x_basis, y_basis, z_basis], grid_dtype=np.float64)
x = domain.grid(0)[:, 0, 0]
y = domain.grid(1)[0, :, 0]
z = domain.grid(2)[0, 0]
field1 = domain.new_field(name='field1')
field2 = domain.new_field(name='field2')

for tind, t in enumerate(ts):
    field1['c'][1, 1, :] = field1_data[tind]
    field2['c'][1, 1, :] = field2_data[tind]
    field1_out = field1['g']
    field2_out = field2['g']

    plt.subplot(2, 2, 1)
    plt.pcolormesh(x, z, field1_out[:, 0, :].T)
    plt.colorbar()
    plt.title(r'{}, $t={:4.2f}$'.format(field1_name, t))
    # plt.xlabel(r'$x$')
    plt.ylabel(r'$z$')

    plt.subplot(2, 2, 2)
    plt.pcolormesh(x, z, field2_out[:, 0, :].T)
    plt.colorbar()
    plt.title(r'{}'.format(field2_name))
    # plt.xlabel(r'$x$')
    # plt.ylabel(r'$z$')

    plt.subplot(2, 2, 3)
    plt.pcolormesh(x, z, field1_out[:, 7, :].T)
    plt.colorbar()
    # plt.title(r'$\Omega_y$')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$z$')

    plt.subplot(2, 2, 4)
    plt.pcolormesh(x, z, field2_out[:, 7, :].T)
    plt.colorbar()
    # plt.title(r'$J_y$')
    plt.xlabel(r'$x$')
    # plt.ylabel(r'$z$')
    # plt.tight_layout()
    plt.savefig(str(plot_dir) + '/snap_{:05d}'.format(tind))
    plt.close()
    print(tind, "/", len(ts)-1)
# plt.show()
