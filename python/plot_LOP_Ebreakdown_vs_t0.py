import numpy as np
from matplotlib import pyplot as plt
import h5py
from configparser import ConfigParser
from pathlib import Path
from dedalus import public as de
import logging
logger = logging.getLogger(__name__)

run = "run_kx-log-ky-lin_MA0.25_Re50_kxmin1e-3_T3tA_Nz256-k300-saveLOPs"
# run = "run_kx-lin-ky-lin_MA1_Re50_kxmin0.1_T20_Nz256-k300-saveLOPs"
kxi = None
kyi = None

if kxi is None and kyi is None:
    get_max_ks = True
else:
    get_max_ks = False

# The .cfg file for this run
config_file = Path("LOPs_3D/config_files/{}.cfg".format(run))
# data dir that was saved to
basedir = Path("LOPs_3D") / config_file.stem
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

gains_file = "LOPs_3D/{}/gains.h5".format(run)
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
    kxs_opt = kxs[kxind]
    kys_opt = kys[kyind]

if get_max_ks:
    for LOP_t0_ind in range(len(LOP_ts)):
        kxky_dir = kxky_dirs[LOP_t0_ind]
        with h5py.File(str(kxky_dir / 'gains_mats.h5'), mode='r') as file:
            # Gs_local = np.array(file['gains/Gs'])
            if LOP_t0_ind == 0:
                LOPs = np.zeros_like(np.array(file['gains/LOPs']))
            LOPs[:, LOP_t0_ind] = np.array(file['gains/LOPs'])[:, LOP_t0_ind]
else:
    kxky_dir = basedir / 'kxi_{:04d}/kyi_{:04d}/'.format(kxi, kyi)
    with h5py.File(str(kxky_dir / 'gains_mats.h5'), mode='r') as file:
        # Gs_local = np.array(file['gains/Gs'])
        LOPs = np.array(file['gains/LOPs'])
    kx = kxs[kxi]
    ky = kys[kyi]

z_basis = de.Chebyshev('z', Nz, interval=(-0.5 * Lz, 0.5 * Lz))
domain = de.Domain([z_basis], grid_dtype=np.complex128)
z = domain.grid(0)
probvars = ['u', 'v', 'w', 'p', 'Ax', 'Ay', 'Az', 'Phi', 'omega_x', 'omega_y', 'Bx', 'By']
fields = [de.Field(domain, name=probvar) for probvar in probvars]
Bz = de.Field(domain, name='Bz')
generic_field = de.Field(domain, name='f')
int_of_field = z_basis.Integrate(generic_field)
energy_components = np.zeros((6, len(LOP_ts)), dtype=np.float64)

for t0_ind in range(len(LOP_ts)):
    if get_max_ks:
        kx = kxs[kxind[t0_ind]]
        ky = kys[kyind[t0_ind]]
    for i, var in enumerate(probvars):
        fields[i]['c'] = LOPs[:, t0_ind][int(i*len(LOPs[:, t0_ind])/len(probvars)):int((i+1)*len(LOPs[:, t0_ind])/len(probvars))]
    for out_i, field_i in enumerate([0, 1, 2, -2, -1]):  # indices in meansquares vs probvars corresponding to u_i and b_i
        if out_i < 3:
            generic_field['g'] = np.abs(fields[field_i]['g'])**2.0
        else:  # for bx, by, bz, make sure to divide by MA^2
            generic_field['g'] = np.abs(fields[field_i]['g']) ** 2.0 / MA2
        energy_components[out_i, t0_ind] = int_of_field.evaluate()['g'][0]
    Bz['c'] = 1j*kx*fields[5]['c'] - 1j*ky*fields[4]['c']  # Bz = dx(Ay) - dy(Ax)
    generic_field['g'] = np.abs(Bz['g'])**2.0 / MA2
    energy_components[-1, t0_ind] = int_of_field.evaluate()['g'][0]

energy_cumsum = np.cumsum(energy_components, axis=0)
colors = ['lightcoral', 'brown', 'red', 'skyblue', 'navy', 'blue']
names = ['u', 'v', 'w', 'bx', 'by', 'bz']
for i in range(6)[::-1]:
    plt.plot(LOP_ts[1:], energy_cumsum[i, 1:], ls='--', c=colors[i], label=names[i])
    if i == 0:
        plt.fill_between(LOP_ts[1:], energy_cumsum[i, 1:], color=colors[i])
    else:
        plt.fill_between(LOP_ts[1:], energy_cumsum[i, 1:], energy_cumsum[i-1, 1:], color=colors[i])
plt.legend()
plt.xlim((0, LOP_ts[-1]))
plt.ylim((0, 1))
plt.xlabel(r'$t_0$')
plt.ylabel(r'Contribution to energy by different fields')
plt.savefig(str(basedir / 'LOP_Ebreakdown_{}.png'.format(run[4:])))
plt.savefig(str(basedir / 'LOP_Ebreakdown_{}.pdf'.format(run[4:])))
