import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import h5py
from configparser import ConfigParser
from pathlib import Path
# import glob

logscale = False
hydro = False
run = "run_Q-k0test"

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
ky_start = params.getfloat('ky_start')
ky_stop = params.getfloat('ky_stop')
ky_num = params.getint('ky_num')
ky_endpoint = params.getboolean('ky_endpoint')

t_start = params.getfloat('t_start')
t_stop = params.getfloat('t_stop')
t_num = params.getint('t_num')
t_endpoint = params.getboolean('t_endpoint')
times = np.linspace(t_start, t_stop, t_num, endpoint=t_endpoint)

k = params.getint('k')
dense = params.getboolean('dense')
Nz = params.getint('Nz')
Lz = params.getfloat('Lz_factor')*np.pi

kys = np.linspace(ky_start, ky_stop, ky_num, endpoint=ky_endpoint)  # the full set of kys to scan over
kxs = np.linspace(kx_start, kx_stop, kx_num, endpoint=kx_endpoint)  # the full set of kxs to scan over
kxdirs = [basedir / 'kxi_{:04d}/'.format(kx_i) for kx_i, kx in enumerate(kxs)]

if hydro:
    gains_file = "hydro/LOPs_3D/{}/gains.h5".format(run)
else:
    gains_file = "LOPs_3D/{}/gains.h5".format(run)
with h5py.File(gains_file, mode='r') as file:
    nEVs = np.array(file['EVs/nEVs'])

if hydro:
    plotname = 'hydro/LOPs_3D/{}/nEVs_kscan_'.format(run)
else:
    plotname = 'LOPs_3D/{}/nEVs_kscan_'.format(run)
plt.pcolormesh(kxs, kys, nEVs)
plt.colorbar()
plt.xlabel(r'$k_x$')
plt.ylabel(r'$k_y$')

if hydro:
    plotname += 'Re{}_'.format(Reynolds)
else:
    plotname += 'MA{}_Re{}_Pm{}_'.format(MA, Reynolds, Pm)
if dense:
    plotname += 'dense'
else:
    plotname += 'k{}'.format(k)
plt.savefig(plotname+'.png')
plt.savefig(plotname+'.pdf')
print("done")
