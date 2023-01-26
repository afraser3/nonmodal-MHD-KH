import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import h5py
from configparser import ConfigParser
from pathlib import Path
# import glob

logscale = True
hydro = False
plotname = '/Users/adfraser/PycharmProjects/nonmodal-MHD-KH/tex/DFD_2022/gains_kmax_scan.pdf'
scale = 0.6
plt.figure(figsize=(scale*7, scale*5))
# for run in ["run_Q", "run_P", "run_L2", "run_U", "run_T"]:
for run in ["run_Q", "run_P", "run_L4", "run_U2"][::-1]:#, "run_T2"]:
    # run = "run_Q"

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
        # kxs = np.array(file['ks/kx'])
        # kys = np.array(file['ks/ky'])
        Gs = np.array(file['gains/Gs'])
        # times = np.array(file['scales/ts'])

    # plotname = '/Users/adfraser/PycharmProjects/nonmodal-MHD-KH/tex/DFD_2022/gains_kmax_'

    # plt.subplot(1, 2, 1)
    plt.plot(times, np.amax(Gs, axis=(0, 1)), label=r'$M_A = {}$'.format(MA))
plt.legend()
plt.xlabel(r'Target time $T$')
# plt.title(r'$G$ maximized over $k_x$ and $k_y$')
plt.ylabel(r'Energy gain $G(T)$')
if logscale:
    plt.yscale("log")
else:
    plt.ylim(ymin=0)
# plt.xlim((times[0], times[-1]))
plt.xlim((0, 800))
plt.ylim((1, 1e5))

plt.tight_layout()
plt.savefig(plotname)
# plt.show()
print("done")
