import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import h5py
from configparser import ConfigParser
from pathlib import Path
# import glob

plot_max = False
logscale = False
hydro = False
# targets = np.array([1.0, 5.0, 7.0, 20.0, 25.0, 50.0])
# targets = np.array([1.0, 5.0, 10.0, 15.0, 25.0, 50.0])
# targets = np.array([1.0, 5.0, 25.0, 50.0, 100.0, 200.0])
# targets = np.array([1.0, 25.0, 100.0, 200.0, 400.0, 800.0])
targets = np.array([1, 25, 150])
run = "run_U2"  # L3/U2

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


plotname = '/Users/adfraser/PycharmProjects/nonmodal-MHD-KH/tex/DFD_2022/gains_kscan_'
if plot_max:
    plotname += 'tmax_'
    if logscale:
        plt.pcolormesh(kxs, kys, np.amax(Gs, axis=2), norm=colors.LogNorm())
    else:
        plt.pcolormesh(kxs, kys, np.amax(Gs, axis=2))
    plt.colorbar()
    plt.xlabel(r'$k_x$')
    plt.ylabel(r'$k_y$')
    # plt.show()
    # plt.savefig('LOPs_3D/{}/gains_tmax_kscan.png'.format(run))
    # plt.savefig('LOPs_3D/{}/gains_tmax_kscan.pdf'.format(run))

else:
    plotname += 'tscan_'
    scale = 0.55
    plt.figure(figsize=(scale*13, scale*4))
    for i, t in enumerate(targets):
        plt.subplot(1, 3, i+1)
        tind = np.argmin(np.abs(times - t))  # index in times where t occurs
        if logscale:
            plt.pcolormesh(kxs, kys, Gs[:, :, tind], norm=colors.LogNorm())
        else:
            plt.pcolormesh(kxs, kys, Gs[:, :, tind])
        if i == 2:
            plt.colorbar(label=r'$G(T)$')
        else:
            plt.colorbar()
        plt.title(r'$T = {}$'.format(t))
        # if i > 3:
        if True:
            plt.xlabel(r'$k_x$ (streamwise)')
        if i == 0:
            plt.ylabel(r'$k_y$ (spanwise)')
    # plt.show()
    # if logscale:
        # plt.savefig('LOPs_3D/{}/gains_tscan_kscan_{}_log.png'.format(run, run))
        # plt.savefig('LOPs_3D/{}/gains_tscan_kscan_{}_log.pdf'.format(run, run))
        # plt.savefig('LOPs_3D/{}/gains_tscan_kscan_MA{}_Re{}_log.png'.format(run, run))
    # else:
        # plt.savefig('LOPs_3D/{}/gains_tscan_kscan_{}_lin.png'.format(run, run))
        # plt.savefig('LOPs_3D/{}/gains_tscan_kscan_{}_lin.pdf'.format(run, run))
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
# plt.savefig(plotname+'.png')
plt.savefig(plotname+'.pdf')
# plt.show()
print("done")
