import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import h5py
from configparser import ConfigParser
from pathlib import Path
# import glob

logscale = False
hydro = False
run = "run_kx-lin-ky-lin_MA0.5_Re250_kxmin0.1_T100_Nz512-k600-saveLOPsAndEvolved"
# run = "run_kx-lin-ky-lin_MA1_Re250_kxmin0.1_Nz512-k600-saveLOPsAndEvolved"

# Parse the .cfg filename passed to script
if hydro:
    config_file = Path("hydro/LOPs_3D/config_files/{}.cfg".format(run))
else:
    config_file = Path("LOPs_3D/config_files/{}.cfg".format(run))

# Parse settings from .cfg file and save to variables
runconfig = ConfigParser()
runconfig.read(str(config_file))

# data dir that was saved to
if hydro:
    basedir = Path("hydro/LOPs_3D") / config_file.stem
else:
    basedir = Path("LOPs_3D") / config_file.stem
# if not basedir.exists():
    # basedir.mkdir(parents=True)

# Parameters
params = runconfig['parameters']
Reynolds = params.getfloat('Reynolds')
if not hydro:
    MA = params.getfloat('MA')
    MA2 = MA**2.0
    Pm = params.getfloat('Pm')
    mReynolds = Pm*Reynolds
kx_log = params.getboolean('kx_log')
kx_start = params.getfloat('kx_start')
kx_stop = params.getfloat('kx_stop')
kx_num = params.getint('kx_num')
kx_endpoint = params.getboolean('kx_endpoint')
add_kx_zero = params.getboolean('add_kx_zero')
ky_start = params.getfloat('ky_start')
ky_stop = params.getfloat('ky_stop')
ky_num = params.getint('ky_num')
ky_endpoint = params.getboolean('ky_endpoint')

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

kys = np.linspace(ky_start, ky_stop, ky_num, endpoint=ky_endpoint)  # the full set of kys to scan over
if kx_log:
    kxs = np.geomspace(kx_start, kx_stop, kx_num, endpoint=kx_endpoint)  # the full set of kxs to scan over
else:
    kxs = np.linspace(kx_start, kx_stop, kx_num, endpoint=kx_endpoint)  # the full set of kxs to scan over
if add_kx_zero and kxs[0] != 0:
    kx_global = np.insert(kxs, 0, 0.0)
    # TODO: pretty sure I forgot to finish implementing this. I think this will screw up the following line (as well as
    #  the kx plots) for runs with add_kx_zero=True, because in the remainder of this file I'm using kxs rather than
    #  kx_global
kxdirs = [basedir / 'kxi_{:04d}/'.format(kx_i) for kx_i, kx in enumerate(kxs)]

if hydro:
    gains_file = "hydro/LOPs_3D/{}/gains.h5".format(run)
else:
    gains_file = "LOPs_3D/{}/gains.h5".format(run)
with h5py.File(gains_file, mode='r') as file:
    # kxs = np.array(file['ks/kx'])
    # kys = np.array(file['ks/ky'])
    Gs = np.array(file['gains/Gs'])
    # times = np.array(file['scales/ts'])

if hydro:
    plotname = 'hydro/LOPs_3D/{}/gains_kmax_{}'.format(run, run)
else:
    plotname = 'LOPs_3D/{}/gains_kmax_{}'.format(run, run)

if np.isnan(Gs[0,0,0]):
    Gs[0, 0] = 0
plt.subplot(2, 2, 1)
plt.plot(times, np.amax(Gs, axis=(0, 1)))
plt.axvline(LOP_ts[4], c='k')
plt.axvline(LOP_ts[10], c='k')
# plt.axvline(LOP_ts[15], c='k')
plt.xlabel(r'$t_0$')
plt.title(r'$G$ maximized over $k_x$ and $k_y$')
plt.ylabel(r'$G(t_0)$')
plt.ylim(ymin=0)
plt.xlim((times[0], times[-1]))

plt.subplot(2, 2, 2)
plt.semilogy(times, np.amax(Gs, axis=(0, 1)))
plt.xlabel(r'$t_0$')
plt.title(r'$G$ maximized over $k_x$ and $k_y$')
plt.ylabel(r'$G(t_0)$')
plt.xlim((times[0], times[-1]))

plt.subplot(2, 2, 3)
# kyind, kxind = np.array([np.unravel_index(np.argmax(Gs[:, :, i]), Gs[:, :, i].shape) for i in range(len(Gs[0, 0]))]).T
Gs_flat = Gs.reshape((-1, Gs.shape[-1]))  # flatten into a 2D array of shape (Gs.shape[0]*Gs.shape[1], Gs.shape[2])
ind_flat = Gs_flat.argmax(axis=0)  # get the index where Gs_flat is maximized for each t
kyind, kxind = np.unravel_index(ind_flat, Gs.shape[:2])  # unflatten/unravel those indices back into their kx, ky shape
plt.plot(times, kxs[kxind], label=r'$k_{x, opt}$')
plt.plot(times, kys[kyind], label=r'$k_{y, opt}$')
plt.axhline(kxs[0], c='C0', ls='--')
plt.axhline(kxs[-1], c='C0', ls='--')
plt.axhline(kys[0], c='C1', ls='--')
plt.axhline(kys[-1], c='C1', ls='--')
plt.legend()
plt.xlabel(r'$t_0$')
plt.ylabel(r'$k$')
plt.ylim(ymin=0)
plt.xlim((times[0], times[-1]))

plt.subplot(2, 2, 4)
plt.semilogy(times, kxs[kxind], label=r'$k_{x, opt}$')
plt.semilogy(times, kys[kyind], label=r'$k_{y, opt}$')
if kxs[0] > 0:
    plt.axhline(kxs[0], c='C0', ls='--')
else:
    plt.axhline(kxs[1], c='C0', ls='--')
plt.axhline(kxs[-1], c='C0', ls='--')
if kys[0] > 0:
    plt.axhline(kys[0], c='C1', ls='--')
else:
    plt.axhline(kys[1], c='C1', ls='--')
plt.axhline(kys[-1], c='C1', ls='--')
plt.legend()
plt.xlabel(r'$t_0$')
plt.ylabel(r'$k$')
plt.xlim((times[0], times[-1]))

# if hydro:
#     plotname += 'Re{}_'.format(Reynolds)
# else:
#     plotname += 'MA{}_Re{}_Pm{}_'.format(MA, Reynolds, Pm)
# if dense:
#     plotname += 'dense'
# else:
#     plotname += 'k{}'.format(k)
plt.tight_layout()
plt.savefig(plotname+'.png')
plt.savefig(plotname+'.pdf')
print("done")
