import numpy as np
from matplotlib import pyplot as plt
plt.style.use('prl.mplstyle')
# import matplotlib.colors as colors
from matplotlib import gridspec as gridspec
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


def make_LOP_plot(ax, LOP2D, x, z):
    return ax.pcolormesh(x, z, LOP2D.T, rasterized=True)


# run_datas_small = [helper_functions.get_Gs_ks_from_run(run) for run in runs_small]
run_datas = [helper_functions.get_Gs_ks_from_run(run) for run in runs]
t0ind1 = LOPt0ind1*run_datas[MAind1]['save_LOPs_istep']
t0ind2 = LOPt0ind2*run_datas[MAind2]['save_LOPs_istep']

# fig = plt.figure(figsize=(18, 4.5), layout="constrained")
fig = plt.figure(figsize=(18, 7.5))#, layout="constrained")
gs0 = gridspec.GridSpec(1, 3, figure=fig)
gs00 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs0[0])
ax1 = fig.add_subplot(gs00[0, :])
ax2 = fig.add_subplot(gs00[1, :])
ax3 = fig.add_subplot(gs00[2, :])
make_G_plot(ax1, runs, run_labels, log=True, ymax=1e4, tmax=100)
ax1.set_ylabel(r'$G(t_0)$')
# ax1.set_title('Small boxes')
# ax1.legend(ncols=2)
ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left', ncols=2, mode="expand", borderaxespad=0.)
ax1.plot(run_datas[MAind1]['times'][t0ind1], run_datas[MAind1]['Gs_opt'][t0ind1], 'D', c='C1')
ax1.plot(run_datas[MAind2]['times'][t0ind2], run_datas[MAind2]['Gs_opt'][t0ind2], 'X', c='C1')
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
ax3.set_xlabel(r'$t_0$')

#########################################################################

gs01 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[1])
gs010 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs01[0])
# gs010 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs01[0])
ax4 = fig.add_subplot(gs010[0, 0])
# ax5 = fig.add_subplot(gs010[1, 0])
ax5 = fig.add_subplot(gs010[0, 1])
ax6 = fig.add_subplot(gs010[1, 0])
ax7 = fig.add_subplot(gs010[1, 1])
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
cmesh4 = make_LOP_plot(ax4, u1[:, 0, :], x16, z16)
cbar4 = fig.colorbar(cmesh4)
cmesh5 = make_LOP_plot(ax5, w1[:, 0, :], x16, z16)
cbar5 = fig.colorbar(cmesh5)
cmesh6 = make_LOP_plot(ax6, Bx1[:, 0, :], x16, z16)
cbar6 = fig.colorbar(cmesh6)
cmesh7 = make_LOP_plot(ax7, Bz1[:, 0, :], x16, z16)
cbar7 = fig.colorbar(cmesh7)

gs011 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs01[1])
# gs011 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs01[1])
ax8 = fig.add_subplot(gs011[0, 0])
# ax9 = fig.add_subplot(gs011[1, 0])
ax9 = fig.add_subplot(gs011[0, 1])
ax10 = fig.add_subplot(gs011[1, 0])
ax11 = fig.add_subplot(gs011[1, 1])
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
cmesh8 = make_LOP_plot(ax8, u2[:, 0, :], x16, z16)
cbar8 = fig.colorbar(cmesh8)
cmesh9 = make_LOP_plot(ax9, w2[:, 0, :], x16, z16)
cbar9 = fig.colorbar(cmesh9)
cmesh10 = make_LOP_plot(ax10, Bx2[:, 0, :], x16, z16)
cbar10 = fig.colorbar(cmesh10)
cmesh11 = make_LOP_plot(ax11, Bz2[:, 0, :], x16, z16)
cbar11 = fig.colorbar(cmesh11)

#########################################################################

gs02 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[2])
gs020 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs02[0])
IVP_data1 = helper_functions.get_IVP_scalar_data(runs[MAind1], LOPt0ind1)
KE1 = IVP_data1['urms']**2 + IVP_data1['vrms']**2 + IVP_data1['wrms']**2
ME1 = (IVP_data1['bxrms']**2 + IVP_data1['byrms']**2 + IVP_data1['bzrms']**2) / run_datas[MAind1]['MA']**2
K_prod1 = IVP_data1['K_prod']
M_prod1 = IVP_data1['M_prod']
ts1 = IVP_data1['scalars_ts']
ax12 = fig.add_subplot(gs020[0])
ax12.plot(ts1, KE1, label=r'$K$')
ax12.plot(ts1, ME1, label=r'$M$')
ax12.set_ylabel(r'Energy')
ax12.set_xlim((0, 100))
ax12.set_ylim(ymin=0)
plt.legend()
ax13 = fig.add_subplot(gs020[1])
ax13.plot(ts1, K_prod1, label=r'$K$ production')
ax13.plot(ts1, M_prod1, label=r'$M$ production')
ax13.set_ylabel(r'Production')
ax13.set_xlabel(r'$t$')
ax13.set_xlim((0, 100))

gs021 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs02[1])
IVP_data2 = helper_functions.get_IVP_scalar_data(runs[MAind2], LOPt0ind2)
KE2 = IVP_data2['urms']**2 + IVP_data2['vrms']**2 + IVP_data2['wrms']**2
ME2 = (IVP_data2['bxrms']**2 + IVP_data2['byrms']**2 + IVP_data2['bzrms']**2) / run_datas[MAind2]['MA']**2
K_prod2 = IVP_data2['K_prod']
M_prod2 = IVP_data2['M_prod']
ts2 = IVP_data2['scalars_ts']
ax14 = fig.add_subplot(gs021[0])
ax14.plot(ts2, KE2, label=r'$K$')
ax14.plot(ts2, ME2, label=r'$M$')
ax14.set_ylabel(r'Energy')
ax14.set_xlim((0, 100))
ax14.set_ylim(ymin=0)
# plt.legend()
ax15 = fig.add_subplot(gs021[1])
ax15.plot(ts2, K_prod2, label=r'$K$ production')
ax15.plot(ts2, M_prod2, label=r'$M$ production')
ax15.set_ylabel(r'Production')
ax15.set_xlabel(r'$t$')
ax15.set_xlim((0, 100))

#########################################################################

for ax in [ax4, ax5, ax8, ax9]:
    ax.tick_params(labelbottom=False)
for ax in [ax5, ax7, ax9, ax11]:
    ax.tick_params(labelleft=False)
for ax in [ax4, ax6, ax8, ax10]:
    ax.set_ylabel(r'$z$')
# for ax in [ax6, ax7, ax10, ax11]:
for ax in [ax10, ax11]:
    ax.set_xlabel(r'$x$')
for ax in [ax4, ax8]:
    ax.set_title(r'$u$')
for ax in [ax5, ax9]:
    ax.set_title(r'$w$')
for ax in [ax6, ax10]:
    ax.set_title(r'$b_x$')
for ax in [ax7, ax11]:
    ax.set_title(r'$b_z$')
ax_topmid = fig.add_subplot(gs010[:])
ax_topmid.axis('off')
ax_topmid.set_title(r'$t_0 = {:2.0f}$ LOP'.format(run_datas[MAind1]['times'][t0ind1]))
ax_botmid = fig.add_subplot(gs011[:])
ax_botmid.axis('off')
ax_botmid.set_title(r'$t_0 = {:2.0f}$ LOP'.format(run_datas[MAind2]['times'][t0ind2]))

plt.suptitle('Small boxes', fontsize='xx-large')

plt.tight_layout()
# fname = '../../tex/manuscript/figures/smallbox_megaplot_kymin5'
fname = '../../tex/manuscript/figures/smallbox_megaplot4'
# fname = '../../tex/manuscript/figures/bigbox_megaplot'
plt.savefig(fname+'.pdf')
plt.savefig(fname+'.png')
plt.savefig(fname+'.eps')
# plt.show()
