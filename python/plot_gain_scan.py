import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import h5py
from configparser import ConfigParser
from pathlib import Path
# import glob

plot_max = False
logscale = True
hydro = False
plot_condition_number = True
# targets = np.array([1.0, 5.0, 7.0, 20.0, 25.0, 50.0])
# targets = np.array([1.0, 5.0, 10.0, 15.0, 25.0, 50.0])
# targets = np.array([1.0, 5.0, 25.0, 50.0, 100.0, 200.0])
# targets = np.array([1.0, 25.0, 100.0, 200.0, 400.0, 800.0])
run = "run_kx-log-ky-lin_MA1_Re50_kxmin1e-8_Lz15pi_T1e4_Nz256-k300-saveLOPs"

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
target_inds = np.geomspace(1, t_num-1, 6).astype(int)

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
    kxs = np.insert(kxs, 0, 0.0)
kxdirs = [basedir / 'kxi_{:04d}/'.format(kx_i) for kx_i, kx in enumerate(kxs)]

if hydro:
    gains_file = "hydro/LOPs_3D/{}/gains.h5".format(run)
else:
    gains_file = "LOPs_3D/{}/gains.h5".format(run)
with h5py.File(gains_file, mode='r') as file:
    # kxs = np.array(file['ks/kx'])
    # kys = np.array(file['ks/ky'])
    Gs = np.array(file['gains/Gs'])
    if plot_condition_number:
        cond_num = np.array(file['EVs/cond'])
    # times = np.array(file['scales/ts'])
print(np.shape(cond_num))

if hydro:
    plotname = 'hydro/LOPs_3D/{}/gains_kscan_'.format(run)
else:
    plotname = 'LOPs_3D/{}/gains_kscan_'.format(run)
if plot_max:
    plotname += 'tmax_'
    if add_kx_zero:
        if logscale:
            plt.pcolormesh(kxs[1:], kys, np.amax(Gs, axis=2)[1:], norm=colors.LogNorm())
        else:
            plt.pcolormesh(kxs[1:], kys, np.amax(Gs, axis=2)[1:])
    else:
        if logscale:
            plt.pcolormesh(kxs, kys, np.amax(Gs, axis=2), norm=colors.LogNorm())
        else:
            plt.pcolormesh(kxs, kys, np.amax(Gs, axis=2))
    plt.colorbar()
    plt.xlabel(r'$k_x$')
    plt.ylabel(r'$k_y$')
    if kx_log:
        plt.xscale('log')
    # plt.show()
    # plt.savefig('LOPs_3D/{}/gains_tmax_kscan.png'.format(run))
    # plt.savefig('LOPs_3D/{}/gains_tmax_kscan.pdf'.format(run))

else:
    plotname += 'tscan_'
    if plot_condition_number:
        plt.figure(figsize=(10, 16))
    else:
        plt.figure(figsize=(10, 12))
    # for i, t in enumerate(targets):
    for i, tind in enumerate(target_inds):
        if plot_condition_number:
            plt.subplot(4, 2, i+1)
        else:
            plt.subplot(3, 2, i + 1)
        # tind = np.argmin(np.abs(times - t))  # index in times where t occurs
        t = times[tind]
        if add_kx_zero:
            if logscale:
                plt.pcolormesh(kxs[1:], kys, Gs[:, 1:, tind], norm=colors.LogNorm())
            else:
                plt.pcolormesh(kxs[1:], kys, Gs[:, 1:, tind])
        else:
            if logscale:
                plt.pcolormesh(kxs, kys, Gs[:, :, tind], norm=colors.LogNorm())
            else:
                plt.pcolormesh(kxs, kys, Gs[:, :, tind])
        plt.colorbar()
        plt.title(r'$t_0 = {}$'.format(t))
        if kx_log:
            plt.xscale('log')
        if i > 3 and not plot_condition_number:
            plt.xlabel(r'$k_x$')
        if i % 2 == 0:
            plt.ylabel(r'$k_y$')
    if plot_condition_number:
        plt.subplot(4, 2, 7)
        if add_kx_zero:
            if logscale:
                plt.pcolormesh(kxs[1:], kys, cond_num[:, 1:], norm=colors.LogNorm())
            else:
                plt.pcolormesh(kxs[1:], kys, cond_num[:, 1:])
        else:
            if logscale:
                plt.pcolormesh(kxs, kys, cond_num[:, :], norm=colors.LogNorm())
            else:
                plt.pcolormesh(kxs, kys, cond_num[:, :])
        plt.colorbar()
        if kx_log:
            plt.xscale('log')
        plt.xlabel(r'$k_x$')
        plt.title('Condition number')

        plt.subplot(4, 2, 8)
        if add_kx_zero:
            if logscale:
                plt.pcolormesh(kxs[1:], kys[2:], cond_num[2:, 1:], norm=colors.LogNorm())
            else:
                plt.pcolormesh(kxs[1:], kys[2:], cond_num[2:, 1:])
        else:
            if logscale:
                plt.pcolormesh(kxs, kys[2:], cond_num[2:, :], norm=colors.LogNorm())
            else:
                plt.pcolormesh(kxs, kys[2:], cond_num[2:, :])
        plt.colorbar()
        if kx_log:
            plt.xscale('log')
        plt.xlabel(r'$k_x$')
        plt.title('Cond num (cut off 2 kys)')
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
plt.savefig(plotname+'.png')
plt.savefig(plotname+'.pdf')
print("done")
