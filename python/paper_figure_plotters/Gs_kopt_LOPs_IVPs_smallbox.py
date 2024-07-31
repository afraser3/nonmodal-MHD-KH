import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec as gridspec
plt.style.use('prl.mplstyle')
import matplotlib as mpl
mpl.rcParams
# import matplotlib.colors as colors
import h5py
from configparser import ConfigParser
from pathlib import Path
from dedalus import public as de
import helper_functions


# run_dict = helper_functions.run_dict_pleiades
run_dict = helper_functions.run_dict
runs = [run_dict[('small', MA)] for MA in [0.5, 1, 1.1, 40]]
run_labels = [r'$M_A = 0.5$',
              r'$M_A = 1$',
              r'$M_A = 1.1$',
              r'$M_A = 40$']
# runs_big = [run_dict[('big', 0.1)], run_dict[('big', 0.5)], run_dict[('big', 1)], run_dict[('big', 1.1, 200)]]
# runs = [run_dict[('big', 0.1)], run_dict[('big', 0.5)], run_dict[('big', 1)], run_dict[('big', 1.1, 200)]]
# run_labels = [r'$M_A = 0.1$',
#               r'$M_A = 0.5$',
#               r'$M_A = 1$',
#               r'$M_A = 1.1$']
MAind1 = 0  # 1  # which index in the above list to highlight in a panel in the plot
LOPt0ind1 = 1  # 4
MAind2 = 0  # 1
LOPt0ind2 = 10  # 43  # 9
# MAind1 = 1
# LOPt0ind1 = 1
# MAind2 = 1
# LOPt0ind2 = 11




def make_G_plot(ax, runs, labels, log=False, TA_rescale=False, ymax=None, tmax=None, colors=None, linestyles=None):
    run_outs = [helper_functions.get_Gs_ks_from_run(run) for run in runs]  # TODO: just pass run_outs, not runs
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
                ax.plot(run_outs[ri]['times'], run_outs[ri]['kx_opt'], label='$k_{x, \mathrm{opt}}$', c='grey', ls='-')
                ax.plot(run_outs[ri]['times'], run_outs[ri]['ky_opt'], label='$k_{y, \mathrm{opt}}$', c='grey', ls='--')
            ax.plot(run_outs[ri]['times'], run_outs[ri]['kx_opt'], c=colors[ri], ls='-')
            ax.plot(run_outs[ri]['times'], run_outs[ri]['ky_opt'], c=colors[ri], ls='--')
    if tmax is None:
        ax.set_xlim((0, max([run_out['times'][-1] for run_out in run_outs])))
    else:
        ax.set_xlim((0, tmax))
    ax.set_ylim(ymin=0)


def make_LOP_plot(ax, LOP2D, x, z, cmap):
    return ax.pcolormesh(x, z, LOP2D.T, rasterized=True, cmap=cmap)


# run_datas_small = [helper_functions.get_Gs_ks_from_run(run) for run in runs_small]
run_datas = [helper_functions.get_Gs_ks_from_run(run) for run in runs]
t0ind1 = LOPt0ind1*run_datas[MAind1]['save_LOPs_istep']
t0ind2 = LOPt0ind2*run_datas[MAind2]['save_LOPs_istep']

# fig = plt.figure(figsize=(18, 4.5), layout="constrained")
fig = plt.figure(figsize=(16, 7.5))#, layout="constrained")
gs0 = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[1/3, 2/3])#, wspace=0.3)
gs00 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs0[0], hspace=0.1)
ax1 = fig.add_subplot(gs00[0, :])
ax2 = fig.add_subplot(gs00[1, :])
ax3 = fig.add_subplot(gs00[2, :])
make_G_plot(ax1, runs, run_labels, log=True, ymax=1e4, tmax=100)
ax1.set_ylabel(r'$G_\mathrm{opt}$')
# ax1.set_title('Small boxes')
# ax1.legend(ncols=2)
# ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left', ncols=2, mode="expand", borderaxespad=0.)
# ax1.legend(fontsize='large')
ax1.plot(run_datas[MAind1]['times'][t0ind1], run_datas[MAind1]['Gs_opt'][t0ind1], 'D', c='C0')
ax1.plot(run_datas[MAind2]['times'][t0ind2], run_datas[MAind2]['Gs_opt'][t0ind2], 'X', c='C0')
evalue_kxind2 = run_datas[2]['kxinds'][-1]
evalue_kxopt2 = run_datas[2]['kx_opt'][-1]
evalue_kyind2 = run_datas[2]['kyinds'][-1]
evalue_kyopt2 = run_datas[2]['ky_opt'][-1]
evalue_kxind3 = run_datas[3]['kxinds'][-1]
evalue_kxopt3 = run_datas[3]['kx_opt'][-1]
evalue_kyind3 = run_datas[3]['kyinds'][-1]
evalue_kyopt3 = run_datas[3]['ky_opt'][-1]
evalue2 = np.real(run_datas[2]['evalues'][evalue_kxind2, evalue_kyind2])
evalue3 = np.real(run_datas[3]['evalues'][evalue_kxind3, evalue_kyind3])
times2 = run_datas[2]['times']
times3 = run_datas[3]['times']
ax1.plot(times2[int(len(times2)/4):], run_datas[2]['Gs_opt'][-1]*1.75*np.exp(2*evalue2*(times2[int(len(times2)/4):] - times2[-1])), ls='--', c='C2')
ax1.plot(times3[int(len(times3)/10):], run_datas[3]['Gs_opt'][-1]*1.75*np.exp(2*evalue3*(times3[int(len(times3)/10):] - times3[-1])), ls='--', c='C3')
make_kopt_plot(ax2, runs, run_labels, just_kx=True, tmax=100)
ax2.set_ylabel(r'$k_{x, \mathrm{opt}}$')
make_kopt_plot(ax3, runs, run_labels, just_ky=True, tmax=100)
ax3.set_ylabel(r'$k_{y, \mathrm{opt}}$')
ax3.set_xlabel(r'$t_0$', labelpad=2)
ax3.legend(fontsize='x-large', loc="center right")
# ax1.get_xaxis().set_ticks([])
# ax2.get_xaxis().set_ticks([])
ax1.tick_params(labelbottom=False)
ax2.tick_params(labelbottom=False)

#########################################################################

gs01 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[1], wspace=0.3)
gs010 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs01[0], hspace=0.275)
gs0100 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs010[0], hspace=0.1, wspace=0.3)
# gs010 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs01[0])
ax4 = fig.add_subplot(gs0100[0, 0])
# ax5 = fig.add_subplot(gs010[1, 0])
ax5 = fig.add_subplot(gs0100[0, 1])
ax6 = fig.add_subplot(gs0100[1, 0])
ax7 = fig.add_subplot(gs0100[1, 1])
LOP1, kx1, ky1, run_data1 = helper_functions.get_LOP_from_run(runs[MAind1], LOPt0ind1)
z_basis = de.Chebyshev('z', run_datas[MAind1]['Nz'], interval=(-0.5 * run_datas[MAind1]['Lz'], 0.5 * run_datas[MAind1]['Lz']))
x_basis = de.Fourier('x', 4, interval=(0.0, 2.0 * np.pi / kx1))
if ky1 > 0:
    y_basis = de.Fourier('x', 4, interval=(0.0, 2.0 * np.pi / ky1))
else:
    y_basis = de.Fourier('x', 4, interval=(0.0, 2.0 * np.pi / kx1))
domain = de.Domain([x_basis, y_basis, z_basis], grid_dtype=np.float64)
x = domain.grid(0)[:, 0, 0]
y = domain.grid(1)[0, :, 0]
z = domain.grid(2)[0, 0, :]
x16 = x_basis.grid(scale=16)
y16 = y_basis.grid(scale=16)
z16 = z_basis.grid(scale=16)
u1 = helper_functions.LOP_to_3D_array(domain, LOP1, 'u', kx1, ky1, scale=16)
w1 = helper_functions.LOP_to_3D_array(domain, LOP1, 'w', kx1, ky1, scale=16)
Bx1 = helper_functions.LOP_to_3D_array(domain, LOP1, 'Bx', kx1, ky1, scale=16)
Bz1 = helper_functions.LOP_to_3D_array(domain, LOP1, 'Bz', kx1, ky1, scale=16)
cmesh4 = make_LOP_plot(ax4, u1[:, 0, :], x16, z16, cmap='RdBu')
cbar4 = fig.colorbar(cmesh4, use_gridspec=True)
cmesh5 = make_LOP_plot(ax5, w1[:, 0, :], x16, z16, cmap='RdBu')
cbar5 = fig.colorbar(cmesh5, use_gridspec=True)
cmesh6 = make_LOP_plot(ax6, Bx1[:, 0, :], x16, z16, cmap='PiYG')
cbar6 = fig.colorbar(cmesh6, use_gridspec=True)
cmesh7 = make_LOP_plot(ax7, Bz1[:, 0, :], x16, z16, cmap='PiYG')
cbar7 = fig.colorbar(cmesh7, use_gridspec=True)

gs0101 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs010[1], hspace=0.1, wspace=0.3)
# gs011 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs01[1])
ax8 = fig.add_subplot(gs0101[0, 0])
# ax9 = fig.add_subplot(gs011[1, 0])
ax9 = fig.add_subplot(gs0101[0, 1])
ax10 = fig.add_subplot(gs0101[1, 0])
ax11 = fig.add_subplot(gs0101[1, 1])
LOP2, kx2, ky2, run_data2 = helper_functions.get_LOP_from_run(runs[MAind2], LOPt0ind2)
z_basis2 = de.Chebyshev('z', run_datas[MAind2]['Nz'], interval=(-0.5 * run_datas[MAind2]['Lz'], 0.5 * run_datas[MAind2]['Lz']))
x_basis2 = de.Fourier('x', 4, interval=(0.0, 2.0 * np.pi / kx2))
if ky2 > 0:
    y_basis2 = de.Fourier('x', 4, interval=(0.0, 2.0 * np.pi / ky2))
else:
    y_basis2 = de.Fourier('x', 4, interval=(0.0, 2.0 * np.pi / kx2))
domain2 = de.Domain([x_basis2, y_basis2, z_basis2], grid_dtype=np.float64)
x = domain2.grid(0)[:, 0, 0]
y = domain2.grid(1)[0, :, 0]
z = domain2.grid(2)[0, 0, :]
x16 = x_basis2.grid(scale=16)
y16 = y_basis2.grid(scale=16)
z16 = z_basis2.grid(scale=16)
u2 = helper_functions.LOP_to_3D_array(domain2, LOP2, 'u', kx2, ky2, scale=16)
w2 = helper_functions.LOP_to_3D_array(domain2, LOP2, 'w', kx2, ky2, scale=16)
Bx2 = helper_functions.LOP_to_3D_array(domain2, LOP2, 'Bx', kx2, ky2, scale=16)
Bz2 = helper_functions.LOP_to_3D_array(domain2, LOP2, 'Bz', kx2, ky2, scale=16)
cmesh8 = make_LOP_plot(ax8, u2[:, 0, :], x16, z16, cmap='RdBu')
cbar8 = fig.colorbar(cmesh8, use_gridspec=True)
cmesh9 = make_LOP_plot(ax9, w2[:, 0, :], x16, z16, cmap='RdBu')
cbar9 = fig.colorbar(cmesh9, use_gridspec=True)
cmesh10 = make_LOP_plot(ax10, Bx2[:, 0, :], x16, z16, cmap='PiYG')
cbar10 = fig.colorbar(cmesh10, use_gridspec=True)
cmesh11 = make_LOP_plot(ax11, Bz2[:, 0, :], x16, z16, cmap='PiYG')
cbar11 = fig.colorbar(cmesh11, use_gridspec=True)

#########################################################################

gs011 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs01[1])
gs0110 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs011[0])
IVP_data1 = helper_functions.get_IVP_scalar_data(runs[MAind1], LOPt0ind1)
KE1 = IVP_data1['urms']**2 + IVP_data1['vrms']**2 + IVP_data1['wrms']**2
ME1 = (IVP_data1['bxrms']**2 + IVP_data1['byrms']**2 + IVP_data1['bzrms']**2) / run_datas[MAind1]['MA']**2
K_prod1 = IVP_data1['K_prod']
M_prod1 = IVP_data1['M_prod']
ts1 = IVP_data1['scalars_ts']
ax12 = fig.add_subplot(gs0110[0])
ax12.plot(ts1, KE1, label=r'$K$')
ax12.plot(ts1, ME1, label=r'$M$')
ax12.set_ylabel(r'Energy')
# ax12.set_xlim((0, 100))
ax12.set_xlim((0, 50))
ax12.set_ylim(ymin=0)
ax12.legend(fontsize='x-large', loc="center right")
ax13 = fig.add_subplot(gs0110[1])
ax13.plot(ts1, K_prod1, label=r'$K$ production')
ax13.plot(ts1, M_prod1, label=r'$M$ production')
ax13.set_ylabel(r'Production')
# ax13.set_xlabel(r'$t$')
# ax13.set_xlim((0, 100))
ax13.set_xlim((0, 50))
ax12.tick_params(labelbottom=False)
ax13.tick_params(labelbottom=False)

gs0111 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs011[1])
IVP_data2 = helper_functions.get_IVP_scalar_data(runs[MAind2], LOPt0ind2)
KE2 = IVP_data2['urms']**2 + IVP_data2['vrms']**2 + IVP_data2['wrms']**2
ME2 = (IVP_data2['bxrms']**2 + IVP_data2['byrms']**2 + IVP_data2['bzrms']**2) / run_datas[MAind2]['MA']**2
K_prod2 = IVP_data2['K_prod']
M_prod2 = IVP_data2['M_prod']
ts2 = IVP_data2['scalars_ts']
ax14 = fig.add_subplot(gs0111[0])
ax14.plot(ts2, KE2, label=r'$K$')
ax14.plot(ts2, ME2, label=r'$M$')
ax14.set_ylabel(r'Energy')
# ax14.set_xlim((0, 100))
ax14.set_xlim((0, 50))
ax14.set_ylim(ymin=0)
# plt.legend()
ax15 = fig.add_subplot(gs0111[1])
ax15.plot(ts2, K_prod2, label=r'$K$ production')
ax15.plot(ts2, M_prod2, label=r'$M$ production')
ax15.set_ylabel(r'Production')
ax15.set_xlabel(r'$t$', labelpad=2)
# ax15.set_xlim((0, 100))
ax15.set_xlim((0, 50))
ax14.tick_params(labelbottom=False)

#########################################################################

for ax in [ax4, ax5, ax8, ax9]:
    ax.tick_params(labelbottom=False)
for ax in [ax5, ax7, ax9, ax11]:
    ax.tick_params(labelleft=False)
for ax in [ax4, ax6, ax8, ax10]:
    ax.set_ylabel(r'$z$', labelpad=2)
# for ax in [ax6, ax7, ax10, ax11]:
for ax in [ax10, ax11]:
    ax.set_xlabel(r'$x$', labelpad=2)
for ax in [ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11]:
    ax.set_ylim((-7.5, 7.5))
# ytext = 10
ytext = 5
for ax in [ax4, ax8]:
    # ax.set_title(r'$u$')
    ax.text(2.5, ytext, r'$u$', fontsize='xx-large')
for ax in [ax5, ax9]:
    # ax.set_title(r'$w$')
    ax.text(2.5, ytext, r'$w$', fontsize='xx-large')
for ax in [ax6, ax10]:
    # ax.set_title(r'$b_x$')
    ax.text(2.5, ytext, r'$b_x$', fontsize='xx-large')
for ax in [ax7, ax11]:
    # ax.set_title(r'$b_z$')
    ax.text(2.5, ytext, r'$b_z$', fontsize='xx-large')
ax_topmid = fig.add_subplot(gs0100[:])
ax_topmid.axis('off')
ax_topmid.set_title(r'$t_0 = {:2.0f}$ LOP'.format(run_datas[MAind1]['times'][t0ind1]), pad=10)
ax_botmid = fig.add_subplot(gs0101[:])
ax_botmid.axis('off')
ax_botmid.set_title(r'$t_0 = {:2.0f}$ LOP'.format(run_datas[MAind2]['times'][t0ind2]), pad=10)

# sublabels = [r'$\mathrm{(a)}$', r'$\mathrm{(b)}$', r'$\mathrm{(c)}$', r'$\mathrm{(d)}$']
sublabels = [r'$\mathrm{{({})}}$'.format(i) for i in 'abcdefghijklmno']
# for axi, ax in enumerate(fig.get_axes()):  # this doesn't play nicely with colorbars or my ax_topmid, ax_botmid thing
for axi, ax in enumerate([ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12, ax13, ax14, ax15]):
    # ax.text(-0.2, 1.07, sublabels[axi], horizontalalignment='center', verticalalignment='center',
    #             transform=ax.transAxes, fontsize='xx-large')
    ax.text(1, 0.99, sublabels[axi], horizontalalignment='right', verticalalignment='top',
            transform=ax.transAxes, fontsize='xx-large')#, bbox=dict(boxstyle="square,pad=0.", ec="none", fc="w"))
    ax.tick_params(direction='in')

# plt.suptitle('Small boxes', fontsize='xx-large')

plt.tight_layout()
# fname = '../../tex/manuscript/figures/smallbox_megaplot_kymin5'
fname = '../../tex/manuscript/figures/smallbox_megaplot'
# fname = '../../tex/manuscript/figures/bigbox_megaplot'
plt.savefig(fname+'.pdf')
plt.savefig(fname+'.png')
# plt.savefig(fname+'.eps')
plt.show()
