import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import h5py
from configparser import ConfigParser
from pathlib import Path
# import glob

logscale = False
run1 = "run_kx-log-ky-lin_MA1_Re50_kxmin1e-3_T1e4_Nz256-k300-saveLOPs"
run2 = "run_kx-log-ky-lin_MA1_Re50_kxmin1e-4_T1e4_Nz256-k300-saveLOPs"
basedir1 = Path("LOPs_3D") / run1

label1 = r'$k_{x, min} = 10^{-3}$'
label2 = r'$k_{x, min} = 10^{-4}$'
tmax = 500
ymax = 3000

def get_Gs_ks_from_run(run):
    run_data = {}  # take everything we think we might want and throw it into a dict

    # Parse the .cfg filename passed to script
    config_file = Path("LOPs_3D/config_files/{}.cfg".format(run))

    # Parse settings from .cfg file and save to variables
    runconfig = ConfigParser()
    runconfig.read(str(config_file))

    # data dir that was saved to
    basedir = Path("LOPs_3D") / config_file.stem
    if not basedir.exists():
        basedir.mkdir(parents=True)

    # Parameters
    params = runconfig['parameters']
    Reynolds = params.getfloat('Reynolds')
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

    run_data['times'] = times

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

    run_data['kys'] = kys
    run_data['kxs'] = kxs
    kxdirs = [basedir / 'kxi_{:04d}/'.format(kx_i) for kx_i, kx in enumerate(kxs)]

    gains_file = "LOPs_3D/{}/gains.h5".format(run)
    with h5py.File(gains_file, mode='r') as file:
        # kxs = np.array(file['ks/kx'])
        # kys = np.array(file['ks/ky'])
        Gs = np.array(file['gains/Gs'])
        # times = np.array(file['scales/ts'])
    if np.isnan(Gs[0, 0, 0]):
        Gs[0, 0] = 0
    run_data['Gs_global'] = Gs
    run_data['Gs_opt'] = np.amax(Gs, axis=(0, 1))
    Gs_flat = Gs.reshape((-1, Gs.shape[-1]))  # flatten into a 2D array of shape (Gs.shape[0]*Gs.shape[1], Gs.shape[2])
    ind_flat = Gs_flat.argmax(axis=0)  # get the index where Gs_flat is maximized for each t
    kyind, kxind = np.unravel_index(ind_flat, Gs.shape[:2])  # unflatten/unravel those indices back into their kx, ky shape
    run_data['kx_opt'] = kxs[kxind]
    run_data['ky_opt'] = kys[kyind]
    return run_data

plotname = 'LOPs_3D/{}/compare_gains_with_{}'.format(run1, run2)
run_out1 = get_Gs_ks_from_run(run1)
run_out2 = get_Gs_ks_from_run(run2)

plt.subplot(2, 2, 1)
plt.plot(run_out1['times'], run_out1['Gs_opt'], label=label1, ls='-', c='C0')
plt.plot(run_out2['times'], run_out2['Gs_opt'], label=label2, ls='--', c='C1')
plt.xlabel(r'$t_0$')
plt.title(r'$G$ maximized over $k_x$ and $k_y$')
plt.ylabel(r'$G(t_0)$')
if ymax is None:
    plt.ylim(ymin=0)
else:
    plt.ylim((0, ymax))
if tmax is None:
    plt.xlim((0, max(run_out1['times'][-1], run_out2['times'][-1])))
else:
    plt.xlim((0, tmax))
plt.legend()

plt.subplot(2, 2, 3)
plt.semilogy(run_out1['times'], run_out1['Gs_opt'], label=label1, ls='-', c='C0')
plt.semilogy(run_out2['times'], run_out2['Gs_opt'], label=label2, ls='--', c='C1')
plt.xlabel(r'$t_0$')
# plt.title(r'$G$ maximized over $k_x$ and $k_y$')
plt.ylabel(r'$G(t_0)$')
if tmax is None:
    plt.xlim((0, max(run_out1['times'][-1], run_out2['times'][-1])))
else:
    plt.xlim((0, tmax))
if ymax is None:
    plt.ylim(ymin=1)
else:
    plt.ylim((1, ymax))

plt.subplot(2, 2, 2)
plt.plot(run_out1['times'], run_out1['Gs_opt'], label=label1, ls='-', c='C0')
plt.plot(run_out2['times'], run_out2['Gs_opt'], label=label2, ls='--', c='C1')
plt.xlabel(r'$t_0$')
plt.title(r'$G$ maximized over $k_x$ and $k_y$')
plt.ylabel(r'$G(t_0)$')
plt.ylim(ymin=0)
plt.xlim((0, max(run_out1['times'][-1], run_out2['times'][-1])))
plt.legend()

plt.subplot(2, 2, 4)
plt.semilogy(run_out1['times'], run_out1['Gs_opt'], label=label1, ls='-', c='C0')
plt.semilogy(run_out2['times'], run_out2['Gs_opt'], label=label2, ls='--', c='C1')
plt.xlabel(r'$t_0$')
# plt.title(r'$G$ maximized over $k_x$ and $k_y$')
plt.ylabel(r'$G(t_0)$')
plt.ylim(ymin=1)
plt.xlim((0, max(run_out1['times'][-1], run_out2['times'][-1])))

# plt.subplot(2, 2, 3)
# Gs_flat = Gs.reshape((-1, Gs.shape[-1]))  # flatten into a 2D array of shape (Gs.shape[0]*Gs.shape[1], Gs.shape[2])
# ind_flat = Gs_flat.argmax(axis=0)  # get the index where Gs_flat is maximized for each t
# kyind, kxind = np.unravel_index(ind_flat, Gs.shape[:2])  # unflatten/unravel those indices back into their kx, ky shape
# plt.plot(times, kxs[kxind], label=r'$k_{x, max}$')
# plt.plot(times, kys[kyind], label=r'$k_{y, max}$')
# plt.axhline(kxs[0], c='C0', ls='--')
# plt.axhline(kxs[-1], c='C0', ls='--')
# plt.axhline(kys[0], c='C1', ls='--')
# plt.axhline(kys[-1], c='C1', ls='--')
# plt.legend()
# plt.xlabel(r'$t_0$')
# plt.ylabel(r'$k$')
# plt.ylim(ymin=0)
# plt.xlim((times[0], times[-1]))
#
# plt.subplot(2, 2, 4)
# plt.semilogy(times, kxs[kxind], label=r'$k_{x, max}$')
# plt.semilogy(times, kys[kyind], label=r'$k_{y, max}$')
# if kxs[0] > 0:
#     plt.axhline(kxs[0], c='C0', ls='--')
# else:
#     plt.axhline(kxs[1], c='C0', ls='--')
# plt.axhline(kxs[-1], c='C0', ls='--')
# if kys[0] > 0:
#     plt.axhline(kys[0], c='C1', ls='--')
# else:
#     plt.axhline(kys[1], c='C1', ls='--')
# plt.axhline(kys[-1], c='C1', ls='--')
# plt.legend()
# plt.xlabel(r'$t_0$')
# plt.ylabel(r'$k$')
# plt.xlim((times[0], times[-1]))
#
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
