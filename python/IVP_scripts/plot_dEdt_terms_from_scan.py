import numpy as np
from matplotlib import pyplot as plt
import h5py
from configparser import ConfigParser
from pathlib import Path

run = "run_kx-lin-ky-lin_MA0.5_Re250_kxmin0.1_T100_Nz512-k600-saveLOPsAndEvolved"
# run = "run_kx-log-ky-lin_MA1.01_Re250_kxmin1e-3_T1e4_Nz512-k600-saveLOPsAndEvolved"
LOP_t0_ind = 10
kxi = None
kyi = None

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

with h5py.File(str(out_dir) + '/scalars/scalars_s1/scalars_s1_p0.h5', mode='r') as file:
    urms = np.array(file['tasks/u_rms'])
    vrms = np.array(file['tasks/v_rms'])
    wrms = np.array(file['tasks/w_rms'])
    ts = np.array(file['scales/sim_time'])
    bxrms = np.array(file['tasks/Bx_rms'])
    byrms = np.array(file['tasks/By_rms'])
    bzrms = np.array(file['tasks/Bz_rms'])
    K_prod = np.array(file['tasks/K_shear_prod'])
    M_prod = np.array(file['tasks/M_shear_prod'])
    EM_work = np.array(file['tasks/EM_work'])

KE = np.real(urms[:, 0])**2.0 + np.real(vrms[:, 0])**2.0 + np.real(wrms[:, 0])**2.0
ME = (np.real(bxrms[:, 0])**2.0 + np.real(byrms[:, 0])**2.0 + np.real(bzrms[:, 0])**2.0)/MA2
delta_KE = np.diff(KE)/ts[1]
delta_ME = np.diff(ME)/ts[1]
plt.subplot(3, 1, 1)
plt.axvline(LOP_ts[LOP_t0_ind], c='k')
plt.axhline(0, c='k')
plt.plot(ts, KE, label='KE')
# plt.plot(ts[:-1], delta_KE, label='d(KE)/dt')
plt.plot(ts, ME, label='ME')
# plt.plot(ts[:-1], delta_ME, label='d(ME)/dt')
# plt.plot(ts, np.real(K_prod[:, 0]), label='KE shear prod')
# plt.plot(ts, np.real(M_prod[:, 0]), label='ME shear prod')
# plt.plot(ts, np.real(EM_work[:, 0]), label='EM work')
plt.legend()
plt.xlim((ts[0], ts[-1]))
plt.ylim(ymin=0)

plt.subplot(3, 1, 2)
plt.axvline(LOP_ts[LOP_t0_ind], c='k')
plt.axhline(0, c='k')
# plt.plot(ts, KE, label='KE')
plt.plot(ts[:-1], delta_KE, label='d(KE)/dt')
plt.plot(ts, np.real(K_prod[:, 0]), label='KE shear prod')
plt.plot(ts, np.real(EM_work[:, 0]), label='EM work')
plt.legend(loc='upper right')
plt.xlim((ts[0], ts[-1]))
ymax = np.max(np.abs(np.append(delta_KE, [np.real(K_prod[:, 0]), np.real(EM_work[:, 0])])))
plt.ylim((-ymax, ymax))

plt.subplot(3, 1, 3)
plt.axvline(LOP_ts[LOP_t0_ind], c='k')
plt.axhline(0, c='k')
# plt.plot(ts, ME, label='ME')
plt.plot(ts[:-1], delta_ME, label='d(ME)/dt')
plt.plot(ts, np.real(M_prod[:, 0]), label='ME shear prod')
plt.plot(ts, np.real(-EM_work[:, 0]), label='-(EM work)')
plt.legend(loc='upper right')
plt.xlim((ts[0], ts[-1]))
ymax = np.max(np.abs(np.append(delta_ME, [np.real(M_prod[:, 0]), np.real(-EM_work[:, 0])])))
plt.ylim((-ymax, ymax))
plt.savefig(str(out_dir) + '/energies_prods_{}_t0{:.2e}.png'.format(run, LOP_ts[LOP_t0_ind]))
plt.savefig(str(out_dir) + '/energies_prods_{}_t0{:.2e}.pdf'.format(run, LOP_ts[LOP_t0_ind]))
plt.show()
