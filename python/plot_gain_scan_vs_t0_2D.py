import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import h5py
from configparser import ConfigParser
from pathlib import Path
# import glob

logscale = False
hydro = False
if hydro:
    raise NotImplementedError
run = "run_A"

# Parse the .cfg filename passed to script
if hydro:
    config_file = Path("hydro/LOPs_2D/config_files/{}.cfg".format(run))
else:
    config_file = Path("LOPs_2D/config_files/{}.cfg".format(run))

# Parse settings from .cfg file and save to variables
runconfig = ConfigParser()
runconfig.read(str(config_file))

# data dir that was saved to
if hydro:
    basedir = Path("hydro/LOPs_2D") / config_file.stem
else:
    basedir = Path("LOPs_2D") / config_file.stem
if not basedir.exists():
    basedir.mkdir(parents=True)

# Parameters
params = runconfig['parameters']
Reynolds = params.getfloat('Reynolds')
if not hydro:
    MA = params.getfloat('MA')
    MA2 = MA**2.0
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
times = np.linspace(t_start, t_stop, t_num, endpoint=t_endpoint)

k = params.getint('k')
dense = params.getboolean('dense')
Nz = params.getint('Nz')
Lz = params.getfloat('Lz_factor')*np.pi

kxs = np.linspace(kx_start, kx_stop, kx_num, endpoint=kx_endpoint)  # the full set of kxs to scan over
kxdirs = [basedir / 'kxi_{:04d}/'.format(kx_i) for kx_i, kx in enumerate(kxs)]

if hydro:
    gains_file = "hydro/LOPs_2D/{}/gains.h5".format(run)
else:
    gains_file = "LOPs_2D/{}/gains.h5".format(run)
with h5py.File(gains_file, mode='r') as file:
    # kxs = np.array(file['ks/kx'])
    Gs = np.array(file['gains/Gs'])
    # times = np.array(file['scales/ts'])

if hydro:
    plotname = 'hydro/LOPs_2D/{}/gains_kmax_'.format(run)
else:
    plotname = 'LOPs_2D/{}/gains_kmax_'.format(run)

plt.subplot(1, 2, 1)
plt.plot(times, np.amax(Gs, axis=0))
plt.xlabel(r'$t_0$')
plt.title(r'$G$ maximized over $k_x$')
plt.ylabel(r'$G(t_0)$')
if logscale:
    plt.yscale("log")
else:
    plt.ylim(ymin=0)
plt.xlim((times[0], times[-1]))

plt.subplot(1, 2, 2)
# kyind, kxind = np.array([np.unravel_index(np.argmax(Gs[:, :, i]), Gs[:, :, i].shape) for i in range(len(Gs[0, 0]))]).T
# Gs_flat = Gs.reshape((-1, Gs.shape[-1]))  # flatten into a 2D array of shape (Gs.shape[0]*Gs.shape[1], Gs.shape[2])
# ind_flat = Gs_flat.argmax(axis=0)  # get the index where Gs_flat is maximized for each t
# kyind, kxind = np.unravel_index(ind_flat, Gs.shape[:2])  # unflatten/unravel those indices back into their kx, ky shape
kxind = np.argmax(Gs, axis=0)
plt.plot(times, kxs[kxind], label=r'$k_{x, max}$')
plt.axhline(kxs[0], c='C0', ls='--')
plt.axhline(kxs[-1], c='C0', ls='--')
# plt.legend()
plt.xlabel(r'$t_0$')
plt.ylabel(r'$k_x$')
plt.ylim(ymin=0)
plt.xlim((times[0], times[-1]))

if hydro:
    plotname += 'Re{}_'.format(Reynolds)
else:
    plotname += 'MA{}_Re{}_Pm{}_'.format(MA, Reynolds, Pm)
if dense:
    plotname += 'dense_'
else:
    plotname += 'k{}_'.format(k)
if logscale:
    plotname += 'log'
else:
    plotname += 'lin'
plt.tight_layout()
plt.savefig(plotname+'.png')
plt.savefig(plotname+'.pdf')
print("done")
