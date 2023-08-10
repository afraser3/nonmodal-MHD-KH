"""
For comparing G and k_opt curves, akin to plot_gain_scan_vs_t0.py, but for comparing multiple runs on a single figure
"""

import numpy as np
from matplotlib import pyplot as plt
plt.style.use('prl.mplstyle')
# import matplotlib.colors as colors
from matplotlib import gridspec as gridspec
import h5py
from configparser import ConfigParser
from pathlib import Path

# a dict that maps ('big'/'small', MA) or ('big'/'small', MA, t0_max) tuples to the appropriate run name
run_dict = {('big', 0.1): "run_kx-log-ky-lin_MA0.1_Re250_kxmin1e-3_T3tA_Nz512-k600-saveLOPsAndEvolved",
            ('big', 0.25): "run_kx-log-ky-lin_MA0.25_Re250_kxmin1e-3_T3tA_Nz512-k600-saveLOPsAndEvolved",
            ('big', 0.5): "run_kx-log-ky-lin_MA0.5_Re250_kxmin1e-3_T3tA_Nz512-k600-saveLOPsAndEvolved",
            ('big', 0.9): "run_kx-log-ky-lin_MA0.9_Re250_kxmin1e-3_T3tA_Nz512-k600-saveLOPsAndEvolved",
            ('big', 1): "run_kx-log-ky-lin_MA1_Re250_kxmin1e-3_T3tA_Nz512-k600-saveLOPsAndEvolved",
            ('big', 1.1, 5): "run_kx-log-ky-lin_MA1.1_Re250_kxmin1e-3_T5_Nz512-k600-saveLOPsAndEvolved",
            ('big', 1.1, 200): "run_kx-log-ky-lin_MA1.1_Re250_kxmin1e-3_T200_Nz512-k600-saveLOPsAndEvolved",
            ('big', 40): "run_kx-log-ky-lin_MA40_Re250_kxmin1e-3_T30_Nz512-k600-saveLOPsAndEvolved",
            ('small', 0.1): "run_kx-lin-ky-lin_MA0.1_Re250_kxmin0.1_T100_Nz512-k600-saveLOPsAndEvolved",
            ('small', 0.25): "run_kx-lin-ky-lin_MA0.25_Re250_kxmin0.1_T100_Nz512-k600-saveLOPsAndEvolved",
            ('small', 0.5): "run_kx-lin-ky-lin_MA0.5_Re250_kxmin0.1_T100_Nz512-k600-saveLOPsAndEvolved",
            ('small', 1): "run_kx-lin-ky-lin_MA1_Re250_kxmin0.1_Nz512-k600-saveLOPsAndEvolved",
            ('small', 1.1): "run_kx-lin-ky-lin_MA1.1_Re250_kxmin0.1_T200_Nz512-k600-saveLOPsAndEvolved",
            ('small', 40): "run_kx-lin-ky-lin_MA40_Re250_kxmin0.1_T100_Nz512-k600-saveLOPsAndEvolved"}

# runs = [run_dict[('big', MA)] for MA in [0.1, 0.25, 0.5, 1]]
# labels = [r'$M_A = 0.1$',
#           r'$M_A = 0.25$',
#           r'$M_A = 0.5$',
#           r'$M_A = 1$']
runs_big = [run_dict[('big', MA)] for MA in [0.5, 1, 40]]
runs_big.insert(2, run_dict[('big', 1.1, 200)])
labels_big = [r'$M_A = 0.5$',
              r'$M_A = 1$',
              r'$M_A = 1.1$',
              r'$M_A = 40$']
runs_small = [run_dict[('small', MA)] for MA in [0.5, 1, 1.1, 40]]
labels_small = [r'$M_A = 0.5$',
                r'$M_A = 1$',
                r'$M_A = 1.1$',
                r'$M_A = 40$']
runs_MA1p1 = [run_dict[('small', 1.1)], run_dict[('big', 1.1, 200)]]
labels_MA1p1 = ['small box', 'big box']

tmax = None
ymax = None


def get_Gs_ks_from_run(run):
    run_data = {}  # take everything we think we might want and throw it into a dict

    # Parse the .cfg filename passed to script
    config_file = Path("../LOPs_3D/config_files/{}.cfg".format(run))

    # Parse settings from .cfg file and save to variables
    runconfig = ConfigParser()
    runconfig.read(str(config_file))

    # data dir that was saved to
    basedir = Path("../LOPs_3D") / config_file.stem
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
    run_data['MA'] = MA
    kxdirs = [basedir / 'kxi_{:04d}/'.format(kx_i) for kx_i, kx in enumerate(kxs)]

    gains_file = "../LOPs_3D/{}/gains.h5".format(run)
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


# plotname = 'LOPs_3D/{}/compare_gains_with_{}'.format(run1, run2)
# run_outs_big = [get_Gs_ks_from_run(run) for run in runs_big]
# run_outs_small = [get_Gs_ks_from_run(run) for run in runs_small]


def make_G_plot(ax, runs, labels, log=False, TA_rescale=False, ymax=None, tmax=None, colors=None, linestyles=None):
    run_outs = [get_Gs_ks_from_run(run) for run in runs]  # TODO: just pass run_outs, not runs
    if colors is None:
        colors = ['C{}'.format(ri) for ri in range(len(runs))]
    if linestyles is None:
        linestyles = ['-' for ri in range(len(runs))]
    for ri, run in enumerate(runs):
        if TA_rescale:
            ax.plot(run_outs[ri]['times'] / (2 * np.pi * run_outs[ri]['MA'] / run_outs[ri]['kxs'][0]),
                     run_outs[ri]['Gs_opt'], label=labels[ri], c=colors[ri], ls=linestyles[ri])
        else:
            ax.plot(run_outs[ri]['times'], run_outs[ri]['Gs_opt'], label=labels[ri], c=colors[ri], ls=linestyles[ri])

    if log:
        ax.set_yscale('log')
        ymin = 1
    else:
        ymin = 0
    if tmax is None:
        if TA_rescale:
            ax.set_xlim((0, max([run_out['times'][-1] / (2 * np.pi * run_out['MA'] / run_out['kxs'][0]) for run_out in
                              run_outs])))
        else:
            ax.set_xlim((0, max([run_out['times'][-1] for run_out in run_outs])))

    else:
        ax.set_xlim((0, tmax))
    if ymax is None:
        ax.set_ylim(ymin=ymin)
    else:
        ax.set_ylim((ymin, ymax))


def make_kopt_plot(ax, runs, labels, just_kx=False, just_ky=False, tmax=None, colors=None, linestyles=None):
    run_outs = [get_Gs_ks_from_run(run) for run in runs]
    if colors is None:
        colors = ['C{}'.format(ri) for ri in range(len(runs))]
    if linestyles is None:
        linestyles = ['-' for ri in range(len(runs))]
    for ri, run in enumerate(runs):
        if just_kx:
            ax.plot(run_outs[ri]['times'], run_outs[ri]['kx_opt'], label=labels[ri], c=colors[ri], ls=linestyles[ri])
        elif just_ky:
            ax.plot(run_outs[ri]['times'], run_outs[ri]['ky_opt'], label=labels[ri], c=colors[ri], ls=linestyles[ri])
        else:
            if ri == 0:  # lazy way to make the legend I want
                ax.plot(run_outs[ri]['times'], run_outs[ri]['kx_opt'], label='$k_{x, \mathrm{opt}}$', c='grey', ls='-')
                ax.plot(run_outs[ri]['times'], run_outs[ri]['ky_opt'], label='$k_{y, \mathrm{opt}}$', c='grey', ls='--')
            ax.plot(run_outs[ri]['times'], run_outs[ri]['kx_opt'], c=colors[ri], ls='-')
            ax.plot(run_outs[ri]['times'], run_outs[ri]['ky_opt'], c=colors[ri], ls='--')
    if tmax is None:
        ax.set_xlim((0, max([run_out['times'][-1] for run_out in run_outs])))
    else:
        ax.set_xlim((0, tmax))
    ax.set_ylim(ymin=0)


fig = plt.figure(figsize=(18, 4.5), layout="constrained")
gs0 = gridspec.GridSpec(1, 3, figure=fig)
# gs00 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs0[0])
# ax1 = fig.add_subplot(gs00[:-1, :])
# ax2 = fig.add_subplot(gs00[-1, :-1])
# ax3 = fig.add_subplot(gs00[-1, -1])
gs00 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs0[0])
ax1 = fig.add_subplot(gs00[0, :])
ax2 = fig.add_subplot(gs00[1, :])
ax3 = fig.add_subplot(gs00[2, :])
make_G_plot(ax1, runs_small, labels_small, log=True, ymax=1e4, tmax=100)
ax1.set_ylabel(r'$G(t_0)$')
# ax1.set_xlabel(r'$t_0$')
ax1.set_title('Small boxes')
ax1.legend(ncols=2)
make_kopt_plot(ax2, runs_small, labels_small, just_kx=True, tmax=100)
ax2.set_ylabel(r'$k_{x, \mathrm{opt}}$')
# ax2.set_xlabel(r'$t_0$')
make_kopt_plot(ax3, runs_small, labels_small, just_ky=True, tmax=100)
ax3.set_ylabel(r'$k_{y, \mathrm{opt}}$')
ax3.set_xlabel(r'$t_0$')

# gs01 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs0[1])
# ax4 = fig.add_subplot(gs01[:-1, :])
# ax5 = fig.add_subplot(gs01[-1, :-1])
# ax6 = fig.add_subplot(gs01[-1, -1])
gs01 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs0[1])
ax4 = fig.add_subplot(gs01[0, :])
ax5 = fig.add_subplot(gs01[1, :])
ax6 = fig.add_subplot(gs01[2, :])
make_G_plot(ax4, runs_big, labels_big, log=True, tmax=4e3)
# ax4.set_ylabel(r'$G(t_0)$')
# ax4.set_xlabel(r'$t_0$')
ax4.set_title('Big boxes')
make_kopt_plot(ax5, runs_big, labels_big, just_kx=True, tmax=4e3)
# ax5.set_ylabel(r'$k_{x, \mathrm{opt}}$')
# ax5.set_xlabel(r'$t_0$')
make_kopt_plot(ax6, runs_big, labels_big, just_ky=True, tmax=4e3)
# ax6.set_ylabel(r'$k_{y, \mathrm{opt}}$')
ax6.set_xlabel(r'$t_0$')

# gs02 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs0[2])
# ax7 = fig.add_subplot(gs02[:-1, :])
# ax8 = fig.add_subplot(gs02[-1, :-1])
# ax9 = fig.add_subplot(gs02[-1, -1])
gs02 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs0[2])
ax7 = fig.add_subplot(gs02[0, :])
ax8 = fig.add_subplot(gs02[1, :])
ax9 = fig.add_subplot(gs02[2, :])
make_G_plot(ax7, runs_MA1p1, labels_MA1p1, log=True, colors=['C2', 'C2'], linestyles=['-', '--'])
# ax7.set_ylabel(r'$G(t_0)$')
# ax7.set_xlabel(r'$t_0$')
ax7.legend()
ax7.set_title(r'Both boxes, $M_A = 1.1$')
make_kopt_plot(ax8, runs_MA1p1, labels_MA1p1, just_kx=True, colors=['C2', 'C2'], linestyles=['-', '--'])
# ax8.set_ylabel(r'$k_{x, \mathrm{opt}}$')
# ax8.set_xlabel(r'$t_0$')
make_kopt_plot(ax9, runs_MA1p1, labels_MA1p1, just_ky=True, colors=['C2', 'C2'], linestyles=['-', '--'])
# ax9.set_ylabel(r'$k_{y, \mathrm{opt}}$')
ax9.set_xlabel(r'$t_0$')
# plt.subplot(2, 1, 1)
# make_G_plot(log=True, TA_rescale=False, legend=True)
# plt.xlabel(r'$t_0/T_A$')
# plt.xlabel(r'$t_0$')
# plt.ylabel(r'$G(t_0)$')
# plt.subplot(2, 2, 3)
# make_kopt_plot(just_kx=True, legend=False)
# plt.xlabel(r'$t_0$')
# plt.ylabel(r'$k_{x, \mathrm{opt}}$, $k_{y, \mathrm{opt}}$')
# plt.ylabel(r'$k_{x, \mathrm{opt}}$')
# plt.subplot(2, 2, 4)
# make_kopt_plot(just_ky=True, legend=False)
# plt.xlabel(r'$t_0$')
# plt.ylabel(r'$k_{x, \mathrm{opt}}$, $k_{y, \mathrm{opt}}$')
# plt.ylabel(r'$k_{y, \mathrm{opt}}$')

# plt.tight_layout()
# plt.show()
# plt.savefig('proposed_smallbox_G.png')
fname = '../../tex/manuscript/figures/Gains_megaplot'
plt.savefig(fname+'.pdf')
plt.savefig(fname+'.png')
plt.savefig(fname+'.eps')
