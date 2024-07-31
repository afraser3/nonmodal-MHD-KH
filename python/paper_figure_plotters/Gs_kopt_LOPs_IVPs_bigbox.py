import numpy as np
from matplotlib import pyplot as plt
# import matplotlib.colors as colors
from matplotlib import gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import h5py
from configparser import ConfigParser
from pathlib import Path
# from dedalus import public as de
import helper_functions
plt.style.use('prl.mplstyle')


# run_dict = helper_functions.run_dict_pleiades
run_dict = helper_functions.run_dict
runs = [run_dict[('big', 0.1)], run_dict[('big', 0.25)], run_dict[('big', 0.5)], run_dict[('big', 1)], run_dict[('big', 1.1, 200)]]
run_labels = [r'$M_A = 0.1$',
              r'$M_A = 0.25$',
              r'$M_A = 0.5$',
              r'$M_A = 1$',
              r'$M_A = 1.1$']
colors = ['C4', 'C5', 'C0', 'C1', 'C2']
MAind1 = 1
LOPt0ind1 = 11
MAind2 = 3
LOPt0ind2 = 20




def make_G_plot(ax, runs, labels, log=False, TA_rescale=False, G_rescale=False, ymax=None, tmax=None, colors=None, linestyles=None):
    run_outs = [helper_functions.get_Gs_ks_from_run(run) for run in runs]  # should really just pass run_outs/run_data, not runs
    if colors is None:
        colors = ['C{}'.format(ri) for ri in range(len(runs))]
    if linestyles is None:
        linestyles = ['-' for ri in range(len(runs))]
    for ri, run in enumerate(runs):
        if TA_rescale:
            if G_rescale:
                ax.plot(run_outs[ri]['times'] / (2 * np.pi * run_outs[ri]['MA'] / run_outs[ri]['kxs'][0]),
                        run_outs[ri]['Gs_opt'] / run_outs[ri]['MA'], label=labels[ri], c=colors[ri], ls=linestyles[ri])
            else:
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
    run_outs = [helper_functions.get_Gs_ks_from_run(run) for run in runs]
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
                ax.plot(run_outs[ri]['times'], run_outs[ri]['kx_opt'], label='$k_{x, \mathrm{opt}}$', c='grey', ls='--')
                ax.plot(run_outs[ri]['times'], run_outs[ri]['ky_opt'], label='$k_{y, \mathrm{opt}}$', c='grey', ls='-')
            ax.plot(run_outs[ri]['times'], run_outs[ri]['kx_opt'], c=colors[ri], ls='--')
            ax.plot(run_outs[ri]['times'], run_outs[ri]['ky_opt'], c=colors[ri], ls='-')
    if tmax is None:
        ax.set_xlim((0, max([run_out['times'][-1] for run_out in run_outs])))
    else:
        ax.set_xlim((0, tmax))
    ax.set_ylim(ymin=0)


def make_LOP_plot(ax, LOP2D, x, z):
    return ax.pcolormesh(x, z, LOP2D.T, rasterized=True)


# run_datas_small = [helper_functions.get_Gs_ks_from_run(run) for run in runs_small]
run_datas = [helper_functions.get_Gs_ks_from_run(run) for run in runs]
t0ind1 = LOPt0ind1*run_datas[MAind1]['save_LOPs_istep']
t0ind2 = LOPt0ind2*run_datas[MAind2]['save_LOPs_istep']

# fig = plt.figure(figsize=(18, 4.5), layout="constrained")
fig = plt.figure(figsize=(18, 7.5))#, layout="constrained")
gs0 = gridspec.GridSpec(1, 3, figure=fig)
# gs00 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs0[0])
gs00 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[0])
ax1 = fig.add_subplot(gs00[0, :])
ax2 = fig.add_subplot(gs00[1, :])
# ax3 = fig.add_subplot(gs00[2, :])
make_G_plot(ax1, runs, run_labels, log=True, colors=colors, tmax=7550)  #, tmax=5000)
ax1.set_ylabel(r'$G(t_0)$')
# ax1.set_title('Small boxes')
ax1.legend(ncols=2, fontsize='large')
# ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left', ncols=3, mode="expand", borderaxespad=0.)
ax1.plot(run_datas[MAind1]['times'][t0ind1], run_datas[MAind1]['Gs_opt'][t0ind1], 'D', c='C1')
ax1.plot(run_datas[MAind2]['times'][t0ind2], run_datas[MAind2]['Gs_opt'][t0ind2], 'X', c='C1')
ax1.set_ylim(ymax=3e4)
#ax1ins = inset_axes(ax1, width="50%", height="35%", loc=4, borderpad=1)
#make_G_plot(ax1ins, runs, run_labels, log=False, colors=colors, tmax=2000)
# make_G_plot(ax1ins, runs, run_labels, log=False, colors=colors, TA_rescale=True, G_rescale=False)
# ax1ins.set_xlim((0, 0.2))
# ax1ins.set_ylim((0, 25000))
# ax1ins.set_xlim((1e-3, 0.2))
# ax1ins.set_ylim((1, 25000))
#ax1ins.set_xlim((1e0, 2000))
#ax1ins.set_ylim((1, 25000))
#ax1ins.set_xscale('log')
#ax1ins.set_yscale('log')
ts_compare = np.geomspace(1e0, 1e3)
#ax1ins.loglog(ts_compare, ts_compare**2, '--', c='k')
# make_kopt_plot(ax2, runs, run_labels, just_kx=True, colors=colors, tmax=5000)#, tmax=100)
make_kopt_plot(ax2, runs, run_labels, colors=colors, tmax=7550)  #, tmax=5000)#, tmax=100)
ax2.plot([run_datas[-1]['times'][-1], 7550], [run_datas[-1]['kx_opt'][-1], run_datas[-1]['kx_opt'][-1]], ':', c='C2')
# ax2.axhline(run_datas[-1]['kx_opt'][-1], ls=':', c='C2')
ax2.set_ylabel(r'$k_{x, \mathrm{opt}}, k_{y, \mathrm{opt}}$')
ax2.set_xlabel(r'$t_0$')
# ax2.set_ylim((1e-3, 1))
# ax2.set_yscale('log')
ax2.legend()
# make_kopt_plot(ax3, runs, run_labels, just_ky=True, colors=colors, tmax=5000)#, tmax=100)
# ax3.set_ylabel(r'$k_{y, \mathrm{opt}}$')
# ax3.set_xlabel(r'$t_0$')

#########################################################################

# gs01 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[1])
gs01 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs0[1])
# gs010 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs01[0])
# ax3 = fig.add_subplot(gs010[0, 0])
# ax4 = fig.add_subplot(gs010[1, 0])
ax3 = fig.add_subplot(gs01[0, 0])
ax4 = fig.add_subplot(gs01[1, 0])
LOP_components1, evolved_components1, LOP_ts1 = helper_functions.get_E_breakdown_from_run(runs[MAind1])
LOP_cumsum1 = np.cumsum(LOP_components1, axis=0)
evolved_energy_cumsum1 = np.cumsum(evolved_components1, axis=0)
evolved_energy_cumsum1[:, 1:] = evolved_energy_cumsum1[:, 1:]/evolved_energy_cumsum1[-1, 1:]  # normalize it to 1

cumsum_colors = ['lightcoral', 'brown', 'red', 'skyblue', 'navy', 'blue']
cumsum_names = ['u', 'v', 'w', 'bx', 'by', 'bz']
for i in range(6)[::-1]:
    ax3.plot(LOP_ts1[1:], LOP_cumsum1[i, 1:], ls='-', c=cumsum_colors[i], label=cumsum_names[i])
    if i == 0:
        ax3.fill_between(LOP_ts1[1:], LOP_cumsum1[i, 1:], color=cumsum_colors[i])
    else:
        ax3.fill_between(LOP_ts1[1:], LOP_cumsum1[i, 1:], LOP_cumsum1[i-1, 1:], color=cumsum_colors[i])
# ax3.legend()
ax3.set_xlim((0, LOP_ts1[-1]))
ax3.set_ylim((0, 1))
ax3.axvline(run_datas[MAind1]['times'][t0ind1], c='k')
ax3.tick_params(labelbottom=False, labelleft=False)
# ax3.set_title(r'$t_0 = {:2.0f}$ Energy components'.format(run_datas[MAind1]['times'][t0ind1]))
ax3.set_title(r'{} Energy components'.format(run_labels[MAind1]))
ax3.set_ylabel(r'LOP')
for i in range(6)[::-1]:
    ax4.plot(LOP_ts1[1:], evolved_energy_cumsum1[i, 1:], ls='-', c=cumsum_colors[i], label=cumsum_names[i])
    if i == 0:
        ax4.fill_between(LOP_ts1[1:], evolved_energy_cumsum1[i, 1:], color=cumsum_colors[i])
    else:
        ax4.fill_between(LOP_ts1[1:], evolved_energy_cumsum1[i, 1:], evolved_energy_cumsum1[i-1, 1:], color=cumsum_colors[i])
ax4.set_xlim((0, LOP_ts1[-1]))
ax4.set_ylim((0, 1))
ax4.axvline(run_datas[MAind1]['times'][t0ind1], c='k')
ax4.tick_params(labelleft=False)
ax4.set_ylabel(r'Evolved state')

# gs011 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs01[1])
# gs02 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[2])
gs02 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs0[2])
# ax7 = fig.add_subplot(gs011[0, 0])
# ax8 = fig.add_subplot(gs011[1, 0])
ax6 = fig.add_subplot(gs02[0, 0])
ax7 = fig.add_subplot(gs02[1, 0])
LOP_components2, evolved_components2, LOP_ts2 = helper_functions.get_E_breakdown_from_run(runs[MAind2])
LOP_cumsum2 = np.cumsum(LOP_components2, axis=0)
evolved_energy_cumsum2 = np.cumsum(evolved_components2, axis=0)
evolved_energy_cumsum2[:, 1:] = evolved_energy_cumsum2[:, 1:]/evolved_energy_cumsum2[-1, 1:]  # normalize it to 1

cumsum_colors = ['lightcoral', 'brown', 'red', 'skyblue', 'navy', 'blue']
cumsum_names = ['u', 'v', 'w', 'bx', 'by', 'bz']
for i in range(6)[::-1]:
    ax6.plot(LOP_ts2[1:], LOP_cumsum2[i, 1:], ls='-', c=cumsum_colors[i], label=cumsum_names[i])
    if i == 0:
        ax6.fill_between(LOP_ts2[1:], LOP_cumsum2[i, 1:], color=cumsum_colors[i])
    else:
        ax6.fill_between(LOP_ts2[1:], LOP_cumsum2[i, 1:], LOP_cumsum2[i-1, 1:], color=cumsum_colors[i])
# ax6.legend()
ax6.set_xlim((0, LOP_ts2[-1]))
ax6.set_ylim((0, 1))
ax6.axvline(run_datas[MAind2]['times'][t0ind2], c='k')
ax6.tick_params(labelbottom=False, labelleft=False)
# ax6.set_title(r'$t_0 = {:2.0f}$ Energy components'.format(run_datas[MAind2]['times'][t0ind2]))
ax6.set_title(r'{} Energy components'.format(run_labels[MAind2]))
ax6.set_ylabel(r'LOP')
for i in range(6)[::-1]:
    ax7.plot(LOP_ts2[1:], evolved_energy_cumsum2[i, 1:], ls='-', c=cumsum_colors[i], label=cumsum_names[i])
    if i == 0:
        ax7.fill_between(LOP_ts2[1:], evolved_energy_cumsum2[i, 1:], color=cumsum_colors[i])
    else:
        ax7.fill_between(LOP_ts2[1:], evolved_energy_cumsum2[i, 1:], evolved_energy_cumsum2[i-1, 1:], color=cumsum_colors[i])
ax7.set_xlim((0, LOP_ts2[-1]))
ax7.set_ylim((0, 1))
ax7.axvline(run_datas[MAind2]['times'][t0ind2], c='k')
ax7.tick_params(labelleft=False)
ax7.set_ylabel(r'Evolved state')


#########################################################################

ax5 = fig.add_subplot(gs01[2, 0])
IVP_data1 = helper_functions.get_IVP_scalar_data(runs[MAind1], LOPt0ind1)
KE1 = IVP_data1['urms']**2 + IVP_data1['vrms']**2 + IVP_data1['wrms']**2
ME1 = (IVP_data1['bxrms']**2 + IVP_data1['byrms']**2 + IVP_data1['bzrms']**2) / run_datas[MAind1]['MA']**2
ts1 = IVP_data1['scalars_ts']
# ax5.plot(ts1, KE1+ME1, label=r'$KE+ME$', c='grey', lw=1)
ax5.plot(ts1, KE1+ME1, label=r'$KE+ME$', c='k')
ax5.plot(run_datas[MAind1]['times'], run_datas[MAind1]['Gs_opt'], c='grey', lw=1, label=r'$G(t_0)$')
# ax5.plot(times, Gs_local, c='k', ls='--', label=r'$G(t_0)$ for this kx/ky')
# ax5.plot(run_datas[MAind1]['times'][t0ind1], run_datas[MAind1]['Gs_opt'][t0ind1], 'D', c='C1')
ax5.axvline(run_datas[MAind1]['times'][t0ind1], c='k')
ax5.plot(ts1, KE1, label=r'$KE$')
ax5.plot(ts1, ME1, label=r'$ME$')
ax5.set_yscale('log')
# ax5.set_xlabel(r'$t, t_0$')
# ax5.set_ylabel(r'$G, E$')
# ax5.set_xlabel(r'$t$')
ax5.set_ylabel(r'$E$')
ax5.legend(ncols=2, fontsize='large')
ax5.set_xlim((ts1[0], ts1[-1]))

ax8 = fig.add_subplot(gs02[2, 0])
IVP_data2 = helper_functions.get_IVP_scalar_data(runs[MAind2], LOPt0ind2)
KE2 = IVP_data2['urms']**2 + IVP_data2['vrms']**2 + IVP_data2['wrms']**2
ME2 = (IVP_data2['bxrms']**2 + IVP_data2['byrms']**2 + IVP_data2['bzrms']**2) / run_datas[MAind2]['MA']**2
ts2 = IVP_data2['scalars_ts']
ax8.plot(run_datas[MAind2]['times'], run_datas[MAind2]['Gs_opt'], c='grey', lw=1, label=r'$G(t_0)$')
ax8.plot(ts2, KE2+ME2, label=r'$KE+ME$', c='k')
# ax5.plot(times, Gs_local, c='k', ls='--', label=r'$G(t_0)$ for this kx/ky')
# ax5.plot(run_datas[MAind1]['times'][t0ind1], run_datas[MAind1]['Gs_opt'][t0ind1], 'D', c='C1')
ax8.axvline(run_datas[MAind2]['times'][t0ind2], c='k')
ax8.plot(ts2, KE2, label=r'$KE$')
ax8.plot(ts2, ME2, label=r'$ME$')
ax8.set_yscale('log')
# ax8.set_xlabel(r'$t, t_0$')
# ax8.set_ylabel(r'$G, E$')
# ax8.set_xlabel(r'$t$')
ax8.set_ylabel(r'$E$')
ax8.set_xlim((ts2[0], ts2[-1]))

# ax5.semilogy(scalars_ts, Gs_local[0] * np.exp(scalars_ts * 2.0 * np.max(np.real(eigenvalues[indx]))) * Gs_local[
#     -1] / np.exp(scalars_ts[-1] * 2.0 * np.max(np.real(eigenvalues[indx]))), '--', c='grey',
#              label=r'eigenvalue')
# ax5.legend(ncols=3)

# gs02 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[2])
# gs020 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs02[0])
# IVP_data1 = helper_functions.get_IVP_scalar_data(runs[MAind1], LOPt0ind1)
# KE1 = IVP_data1['urms']**2 + IVP_data1['vrms']**2 + IVP_data1['wrms']**2
# ME1 = (IVP_data1['bxrms']**2 + IVP_data1['byrms']**2 + IVP_data1['bzrms']**2) / run_datas[MAind1]['MA']**2
# K_prod1 = IVP_data1['K_prod']
# M_prod1 = IVP_data1['M_prod']
# ts1 = IVP_data1['scalars_ts']
# ax12 = fig.add_subplot(gs020[0])
# ax12.plot(ts1, KE1, label=r'$K$')
# ax12.plot(ts1, ME1, label=r'$M$')
# ax12.set_ylabel(r'Energy')
# # ax12.set_xlim((0, 100))
# # ax12.set_ylim(ymin=0)
# ax12.set_yscale('log')
# plt.legend()
# ax13 = fig.add_subplot(gs020[1])
# ax13.plot(ts1, K_prod1, label=r'$K$ production')
# ax13.plot(ts1, M_prod1, label=r'$M$ production')
# ax13.set_ylabel(r'Production')
# ax13.set_xlabel(r'$t$')
# # ax13.set_xlim((0, 100))
#
# gs021 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs02[1])
# IVP_data2 = helper_functions.get_IVP_scalar_data(runs[MAind2], LOPt0ind2)
# KE2 = IVP_data2['urms']**2 + IVP_data2['vrms']**2 + IVP_data2['wrms']**2
# ME2 = (IVP_data2['bxrms']**2 + IVP_data2['byrms']**2 + IVP_data2['bzrms']**2) / run_datas[MAind2]['MA']**2
# K_prod2 = IVP_data2['K_prod']
# M_prod2 = IVP_data2['M_prod']
# ts2 = IVP_data2['scalars_ts']
# ax14 = fig.add_subplot(gs021[0])
# ax14.plot(ts2, KE2, label=r'$K$')
# ax14.plot(ts2, ME2, label=r'$M$')
# ax14.set_ylabel(r'Energy')
# # ax14.set_xlim((0, 100))
# # ax14.set_ylim(ymin=0)
# ax14.set_yscale('log')
# # plt.legend()
# ax15 = fig.add_subplot(gs021[1])
# ax15.plot(ts2, K_prod2, label=r'$K$ production')
# ax15.plot(ts2, M_prod2, label=r'$M$ production')
# ax15.set_ylabel(r'Production')
# ax15.set_xlabel(r'$t$')
# # ax15.set_xlim((0, 100))

#########################################################################

# for ax in [ax3, ax4, ax8, ax9]:
#     ax.tick_params(labelbottom=False)
# for ax in [ax4, ax8, ax9, ax11]:
#     ax.tick_params(labelleft=False)
# for ax in [ax3, ax7, ax8, ax10]:
#     ax.set_ylabel(r'$z$')
# for ax in [ax10, ax11]:
#     ax.set_xlabel(r'$x$')
# for ax in [ax3, ax8]:
#     ax.set_title(r'$u$')
# for ax in [ax4, ax9]:
#     ax.set_title(r'$w$')
# for ax in [ax7, ax10]:
#     ax.set_title(r'$b_x$')
# for ax in [ax8, ax11]:
#     ax.set_title(r'$b_z$')
# ax_topmid = fig.add_subplot(gs010[:])
# ax_topmid.axis('off')
# ax_topmid.set_title(r'$t_0 = {:2.0f}$ LOP'.format(run_datas[MAind1]['times'][t0ind1]))
# ax_botmid = fig.add_subplot(gs011[:])
# ax_botmid.axis('off')
# ax_botmid.set_title(r'$t_0 = {:2.0f}$ LOP'.format(run_datas[MAind2]['times'][t0ind2]))
#
# plt.suptitle('Small boxes', fontsize='xx-large')

plt.tight_layout()
# fname = '../../tex/manuscript/figures/smallbox_megaplot_kymin5'
# fname = '../../tex/manuscript/figures/smallbox_megaplot3'
fname = '../../tex/manuscript/figures/bigbox_megaplot2-noinset'
# plt.savefig(fname+'.pdf')
# plt.savefig(fname+'.png')
# plt.savefig(fname+'.eps')
plt.show()
