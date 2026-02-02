import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import gridspec as gridspec
import h5py
from configparser import ConfigParser
from pathlib import Path
from dedalus import public as de
import logging
logger = logging.getLogger(__name__)

# run = "run_kx-lin-ky-lin_MA0.5_Re250_kxmin0.1_T100_Nz512-k600-saveLOPsAndEvolved"
run = "run_kx-log-ky-lin_MA0.25_Re250_kxmin1e-3_T3tA_Nz512-k600-saveLOPsAndEvolved"
LOP_t0_ind = 11
kxi = None
kyi = None

if kxi is None and kyi is None:
    get_max_ks = True
else:
    get_max_ks = False

# The .cfg file for this run
config_file = Path("../LOPs_3D/config_files/{}.cfg".format(run))
# data dir that was saved to
basedir = Path("../LOPs_3D") / config_file.stem
# Parse settings from .cfg file and save to variables
runconfig = ConfigParser()
runconfig.read(str(config_file))

# Parameters
params = runconfig['parameters']
Reynolds = params.getfloat('Reynolds')
MA = params.getfloat('MA')
MA2 = MA**2.0
Pm = params.getfloat('Pm')
mReynolds = Pm*Reynolds

t_start = params.getfloat('t_start')
t_stop = params.getfloat('t_stop')
t_num = params.getint('t_num')
t_endpoint = params.getboolean('t_endpoint')
times = np.linspace(t_start, t_stop, t_num, endpoint=t_endpoint)
save_LOPs_istep = params.getint('save_LOPs_istep')
save_evolved_state = params.getboolean('save_evolved_state')
LOP_ts = times[::save_LOPs_istep]
target_time = LOP_ts[LOP_t0_ind]

k = params.getint('k')
dense = params.getboolean('dense')
Nz = params.getint('Nz')
Lz = params.getfloat('Lz_factor')*np.pi

gains_file = "../LOPs_3D/{}/gains.h5".format(run)
with h5py.File(gains_file, mode='r') as file:
    kx_meshgrid = np.array(file['ks/kx'])
    ky_meshgrid = np.array(file['ks/ky'])
    kxs = kx_meshgrid[0]
    kys = ky_meshgrid[:, 0]
    Gs_global = np.array(file['gains/Gs'])
    Gs = Gs_global[:, :, ::save_LOPs_istep]
    # times = np.array(file['scales/ts'])

if get_max_ks:
    Gs_flat = Gs.reshape((-1, Gs.shape[-1]))  # flatten into a 2D array of shape (Gs.shape[0]*Gs.shape[1], Gs.shape[2])
    ind_flat = Gs_flat.argmax(axis=0)  # get the index where Gs_flat is maximized for each t
    kyind, kxind = np.unravel_index(ind_flat, Gs.shape[:2])  # unflatten/unravel those indices back into their kx, ky shape
    kxky_dirs = [basedir / 'kxi_{:04d}/kyi_{:04d}/'.format(kxind[i], kyind[i]) for i in range(len(Gs[0, 0]))]  # a list where each element points to the directory corresponding to the kx/ky that gives max G at each time
    kxky_dir = kxky_dirs[LOP_t0_ind]
    kxi = kxind[LOP_t0_ind]
    kyi = kyind[LOP_t0_ind]
    kx = kxs[kxi]
    ky = kys[kyi]
else:
    kxky_dir = basedir / 'kxi_{:04d}/kyi_{:04d}/'.format(kxi, kyi)
    kx = kxs[kxi]
    ky = kys[kyi]
if get_max_ks:
    out_dir = basedir / 'IVP_kxi{:04d}_kyi{:04d}_t0i{:03d}'.format(kxi, kyi, LOP_t0_ind)
else:
    out_dir = basedir / 'IVP_specific_kxi{:04d}_kyi{:04d}_t0i{:03d}'.format(kxi, kyi, LOP_t0_ind)

with h5py.File(str(kxky_dir / 'gains_mats.h5'), mode='r') as file:
    Gs_local = np.array(file['gains/Gs'])
    LOP = np.array(file['gains/LOPs'])[:, LOP_t0_ind]
    eigenvalues = np.array(file['setup_mats/eigenvalues'])
    indx = np.argmax(np.real(eigenvalues))

# Get LOP and evolved state from IVP profiles
field_names = ['u', 'v', 'w', 'Bx', 'By', 'Bz']
with h5py.File(str(out_dir) + '/profiles_c/profiles_c_s1/profiles_c_s1_p0.h5', mode='r') as file:
    profiles_ts = np.array(file['scales/sim_time'])
    profiles_ts_target_ind = np.argmin(np.abs(profiles_ts-LOP_ts[LOP_t0_ind]))
    field_data_LOP = dict([(field_names[i], np.array(file['tasks/{}'.format(field_names[i])])[0]) for i in range(len(field_names))])
    field_data_target = dict([(field_names[i], np.array(file['tasks/{}'.format(field_names[i])])[profiles_ts_target_ind]) for i in range(len(field_names))])

# Get IVP scalar data
with h5py.File(str(out_dir / 'scalars/scalars_s1/scalars_s1_p0.h5'), mode='r') as file:
    urms = np.array(file['tasks/u_rms'])
    vrms = np.array(file['tasks/v_rms'])
    wrms = np.array(file['tasks/w_rms'])
    scalars_ts = np.array(file['scales/sim_time'])
    bxrms = np.array(file['tasks/Bx_rms'])
    byrms = np.array(file['tasks/By_rms'])
    bzrms = np.array(file['tasks/Bz_rms'])
    K_prod = np.array(file['tasks/K_shear_prod'])
    M_prod = np.array(file['tasks/M_shear_prod'])
    EM_work = np.array(file['tasks/EM_work'])
# Really there should be a factor of 0.5 on the following
KEs = np.real(urms[:, 0])**2.0 + np.real(vrms[:, 0])**2.0 + np.real(wrms[:, 0])**2.0
MEs = (np.real(bxrms[:, 0])**2.0 + np.real(byrms[:, 0])**2.0 + np.real(bzrms[:, 0])**2.0)/MA2
energies = KEs + MEs
delta_KE = np.diff(KEs)/scalars_ts[1]
delta_ME = np.diff(MEs)/scalars_ts[1]

plotname = "IVP-suite_t0_{:.2e}__{}.pdf".format(LOP_ts[LOP_t0_ind], run)
with PdfPages(str(out_dir) + '/' + plotname) as pdf:
    # plot_gain_scan_vs_t0.py stuff
    if True:
        plt.subplot(2, 2, 1)
        plt.plot(times, np.amax(Gs_global, axis=(0, 1)))
        plt.axhline(Gs_local[::save_LOPs_istep][LOP_t0_ind], c='grey', lw=1.0)
        plt.axvline(LOP_ts[LOP_t0_ind], c='grey', lw=1.0)
        plt.xlabel(r'$t_0$')
        plt.title(r'$G$ maximized over $k_x$ and $k_y$')
        plt.ylabel(r'$G(t_0)$')
        plt.ylim(ymin=0)
        plt.xlim((times[0], times[-1]))

        plt.subplot(2, 2, 2)
        plt.semilogy(times, np.amax(Gs_global, axis=(0, 1)))
        plt.axhline(Gs_local[::save_LOPs_istep][LOP_t0_ind], c='grey', lw=1.0)
        plt.axvline(LOP_ts[LOP_t0_ind], c='grey', lw=1.0)
        plt.xlabel(r'$t_0$')
        plt.title(r'$G$ maximized over $k_x$ and $k_y$')
        plt.ylabel(r'$G(t_0)$')
        plt.xlim((times[0], times[-1]))

        plt.subplot(2, 2, 3)
        Gs_global_flat = Gs_global.reshape((-1, Gs_global.shape[-1]))  # flatten into a 2D array of shape (Gs_global.shape[0]*Gs_global.shape[1], Gs_global.shape[2])
        ind_flat = Gs_global_flat.argmax(axis=0)  # get the index where Gs_global_flat is maximized for each t
        kyind_Gs, kxind_Gs = np.unravel_index(ind_flat,
                                        Gs_global.shape[:2])  # unflatten/unravel those indices back into their kx, ky shape
        plt.plot(times, kxs[kxind_Gs], label=r'$k_{x, max}$')
        plt.plot(times, kys[kyind_Gs], label=r'$k_{y, max}$')
        plt.axvline(LOP_ts[LOP_t0_ind], c='grey', lw=1.0)
        plt.axhline(kxs[0], c='C0', ls='--')
        plt.axhline(kxs[-1], c='C0', ls='--')
        plt.axhline(kys[0], c='C1', ls='--')
        plt.axhline(kys[-1], c='C1', ls='--')
        plt.legend()
        plt.xlabel(r'$t_0$')
        plt.ylabel(r'$k$')
        plt.ylim(ymin=0)
        plt.xlim((times[0], times[-1]))

        plt.subplot(2, 2, 4)
        plt.semilogy(times, kxs[kxind_Gs], label=r'$k_{x, max}$')
        plt.semilogy(times, kys[kyind_Gs], label=r'$k_{y, max}$')
        if kxs[0] > 0:
            plt.axhline(kxs[0], c='C0', ls='--')
        else:
            plt.axhline(kxs[1], c='C0', ls='--')
        plt.axhline(kxs[-1], c='C0', ls='--')
        if kys[0] > 0:
            plt.axhline(kys[0], c='C1', ls='--')
        else:
            plt.axhline(kys[1], c='C1', ls='--')
        plt.axhline(kys[-1], c='C1', ls='--')
        plt.axvline(LOP_ts[LOP_t0_ind], c='grey', lw=1.0)
        plt.legend()
        plt.xlabel(r'$t_0$')
        plt.ylabel(r'$k$')
        plt.xlim((times[0], times[-1]))
        plt.tight_layout()
        pdf.savefig()
        plt.close()

    # plot_LOP_Ebreakdown_vs_t0.py stuff
    if True:
        # data grabbing stuff
        if get_max_ks:
            for LOP_t0_i in range(len(LOP_ts)):
                kdir = kxky_dirs[LOP_t0_i]
                with h5py.File(str(kdir / 'gains_mats.h5'), mode='r') as file:
                    if LOP_t0_i == 0:
                        LOPs = np.zeros_like(np.array(file['gains/LOPs']))
                    LOPs[:, LOP_t0_i] = np.array(file['gains/LOPs'])[:, LOP_t0_i]
                    if save_evolved_state:
                        if LOP_t0_i == 0:
                            evolved_states = np.zeros_like(np.array(file['gains/evolved_states']))
                        evolved_states[:, LOP_t0_i] = np.array(file['gains/evolved_states'])[:, LOP_t0_i]
        else:
            with h5py.File(str(kxky_dir / 'gains_mats.h5'), mode='r') as file:
                LOPs = np.array(file['gains/LOPs'])
                if save_evolved_state:
                    evolved_states = np.array(file['gains/evolved_states'])
        z_basis = de.Chebyshev('z', Nz, interval=(-0.5 * Lz, 0.5 * Lz))
        domain = de.Domain([z_basis], grid_dtype=np.complex128)
        z = domain.grid(0)
        probvars = ['u', 'v', 'w', 'p', 'Ax', 'Ay', 'Az', 'Phi', 'omega_x', 'omega_y', 'Bx', 'By']
        fields = [de.Field(domain, name=probvar) for probvar in probvars]
        Bz = de.Field(domain, name='Bz')
        generic_field = de.Field(domain, name='f')
        int_of_field = z_basis.Integrate(generic_field)
        energy_components = np.zeros((6, len(LOP_ts)), dtype=np.float64)
        evolved_energy_components = np.zeros_like(energy_components)

        for LOP_t0_i in range(len(LOP_ts)):
            if get_max_ks:
                this_kx = kxs[kxind[LOP_t0_i]]
                this_ky = kys[kyind[LOP_t0_i]]
            for i, var in enumerate(probvars):
                fields[i]['c'] = LOPs[:, LOP_t0_i][int(i * len(LOPs[:, LOP_t0_i]) / len(probvars)):int(
                    (i + 1) * len(LOPs[:, LOP_t0_i]) / len(probvars))]
            for out_i, field_i in enumerate(
                    [0, 1, 2, -2, -1]):  # indices in meansquares vs probvars corresponding to u_i and b_i
                if out_i < 3:
                    generic_field['g'] = np.abs(fields[field_i]['g']) ** 2.0
                else:  # for bx, by, bz, make sure to divide by MA^2
                    generic_field['g'] = np.abs(fields[field_i]['g']) ** 2.0 / MA2
                energy_components[out_i, LOP_t0_i] = int_of_field.evaluate()['g'][0]
            Bz['c'] = 1j * this_kx * fields[5]['c'] - 1j * this_ky * fields[4]['c']  # Bz = dx(Ay) - dy(Ax)
            generic_field['g'] = np.abs(Bz['g']) ** 2.0 / MA2
            energy_components[-1, LOP_t0_i] = int_of_field.evaluate()['g'][0]
            if save_evolved_state:
                for i, var in enumerate(probvars):
                    fields[i]['c'] = evolved_states[:, LOP_t0_i][
                                     int(i * len(evolved_states[:, LOP_t0_i]) / len(probvars)):int(
                                         (i + 1) * len(evolved_states[:, LOP_t0_i]) / len(probvars))]
                for out_i, field_i in enumerate(
                        [0, 1, 2, -2, -1]):  # indices in meansquares vs probvars corresponding to u_i and b_i
                    if out_i < 3:
                        generic_field['g'] = np.abs(fields[field_i]['g']) ** 2.0
                    else:  # for bx, by, bz, make sure to divide by MA^2
                        generic_field['g'] = np.abs(fields[field_i]['g']) ** 2.0 / MA2
                    evolved_energy_components[out_i, LOP_t0_i] = int_of_field.evaluate()['g'][0]
                Bz['c'] = 1j * this_kx * fields[5]['c'] - 1j * this_ky * fields[4]['c']  # Bz = dx(Ay) - dy(Ax)
                generic_field['g'] = np.abs(Bz['g']) ** 2.0 / MA2
                evolved_energy_components[-1, LOP_t0_i] = int_of_field.evaluate()['g'][0]

        energy_cumsum = np.cumsum(energy_components, axis=0)
        colors = ['lightcoral', 'brown', 'red', 'skyblue', 'navy', 'blue']
        names = ['u', 'v', 'w', 'bx', 'by', 'bz']

        if save_evolved_state:
            evolved_energy_cumsum = np.cumsum(evolved_energy_components, axis=0)
            evolved_energy_cumsum[:, 1:] = evolved_energy_cumsum[:, 1:] / evolved_energy_cumsum[-1,
                                                                          1:]  # normalize it to 1
            plt.subplot(2, 1, 2)
            for i in range(6)[::-1]:
                plt.plot(LOP_ts[1:], evolved_energy_cumsum[i, 1:], ls='--', c=colors[i], label=names[i])
                if i == 0:
                    plt.fill_between(LOP_ts[1:], evolved_energy_cumsum[i, 1:], color=colors[i])
                else:
                    plt.fill_between(LOP_ts[1:], evolved_energy_cumsum[i, 1:], evolved_energy_cumsum[i - 1, 1:],
                                     color=colors[i])
            plt.legend()
            plt.axvline(LOP_ts[LOP_t0_ind], c='k')
            plt.xlim((0, LOP_ts[-1]))
            plt.ylim((0, 1))
            plt.ylabel(r'Evolved state')
            plt.xlabel(r'$t=t_0$')
            plt.subplot(2, 1, 1)
        for i in range(6)[::-1]:
            plt.plot(LOP_ts[1:], energy_cumsum[i, 1:], ls='--', c=colors[i], label=names[i])
            if i == 0:
                plt.fill_between(LOP_ts[1:], energy_cumsum[i, 1:], color=colors[i])
            else:
                plt.fill_between(LOP_ts[1:], energy_cumsum[i, 1:], energy_cumsum[i - 1, 1:], color=colors[i])
        if not save_evolved_state:
            plt.legend()
        plt.axvline(LOP_ts[LOP_t0_ind], c='k')
        plt.xlim((0, LOP_ts[-1]))
        plt.ylim((0, 1))
        plt.xlabel(r'$t_0$')
        plt.title(r'Contribution to energy by different fields')
        plt.ylabel('LOP')
        plt.tight_layout()
        pdf.savefig()
        plt.close()

    # Energies stuff from the bottom of LOP_IVP_3d_from_scan.py
    if True:
        plt.subplot(2, 1, 1)
        plt.plot(times, np.amax(Gs_global, axis=(0, 1)), c='k', label=r'$G(t_0)$ maxed over kx/ky')
        plt.plot(times, Gs_local, c='k', ls='--', label=r'$G(t_0)$ for this kx/ky')
        plt.plot(LOP_ts[LOP_t0_ind], Gs_local[::save_LOPs_istep][LOP_t0_ind], 'x', c='k')
        plt.plot(scalars_ts, energies, label=r'$KE+ME$')
        plt.plot(scalars_ts, KEs, label=r'$KE$')
        plt.plot(scalars_ts, MEs, label=r'$ME$')
        plt.xlabel(r'$t, t_0$')
        plt.ylabel(r'$G, E$')
        plt.xlim((scalars_ts[0], scalars_ts[-1]))
        plt.ylim(ymin=0)
        # plt.legend()

        plt.subplot(2, 1, 2)
        plt.plot(times, np.amax(Gs_global, axis=(0, 1)), c='k', label=r'$G(t_0)$ maxed over kx/ky')
        plt.plot(times, Gs_local, c='k', ls='--', label=r'$G(t_0)$ for this kx/ky')
        plt.plot(LOP_ts[LOP_t0_ind], Gs_local[::save_LOPs_istep][LOP_t0_ind], 'x', c='k')
        plt.plot(scalars_ts, energies, label=r'$KE+ME$')
        plt.plot(scalars_ts, KEs, label=r'$KE$')
        plt.plot(scalars_ts, MEs, label=r'$ME$')
        plt.yscale('log')
        plt.xlabel(r'$t, t_0$')
        plt.ylabel(r'$G, E$')
        plt.xlim((scalars_ts[0], scalars_ts[-1]))
        plt.semilogy(scalars_ts, Gs_local[0] * np.exp(scalars_ts * 2.0 * np.max(np.real(eigenvalues[indx]))) * Gs_local[
            -1] / np.exp(scalars_ts[-1] * 2.0 * np.max(np.real(eigenvalues[indx]))), '--', c='grey',
                     label=r'eigenvalue')
        plt.legend(ncols=3)
        plt.tight_layout()
        pdf.savefig()
        plt.close()

    # scalars stuff from plot_E_breakdown_from_scan.py
    if True:
        plt.subplot(2, 3, 1)
        plt.plot(scalars_ts, np.real(urms[:, 0]), label='u_rms')
        plt.plot(scalars_ts, np.real(vrms[:, 0]), label='v_rms')
        plt.plot(scalars_ts, np.real(wrms[:, 0]), label='w_rms')
        plt.legend()
        plt.xlim((scalars_ts[0], scalars_ts[-1]))
        # plt.xlim((0, 20))
        plt.ylim(ymin=0)
        plt.title('lin-lin')
        plt.ylabel(r'$\mathbf{u}$ components')

        plt.subplot(2, 3, 4)
        plt.plot(scalars_ts, np.real(bxrms[:, 0]), label='bx_rms')
        plt.plot(scalars_ts, np.real(byrms[:, 0]), label='by_rms')
        plt.plot(scalars_ts, np.real(bzrms[:, 0]), label='bz_rms')
        plt.legend()
        plt.xlim((scalars_ts[0], scalars_ts[-1]))
        # plt.xlim((0, 20))
        plt.ylim(ymin=0)
        plt.ylabel(r'$\mathbf{b}$ components')
        plt.xlabel(r'$t$')

        plt.subplot(2, 3, 2)
        plt.loglog(scalars_ts, np.real(urms[:, 0]), label='u_rms')
        plt.loglog(scalars_ts, np.real(vrms[:, 0]), label='v_rms')
        plt.loglog(scalars_ts, np.real(wrms[:, 0]), label='w_rms')
        # plt.xlim(xmax=20)
        plt.loglog([scalars_ts[10], scalars_ts[-1]], [scalars_ts[10] * urms[10, 0], scalars_ts[-1] * urms[10, 0]], '--', c='k')
        # plt.loglog([scalars_ts[1], scalars_ts[-1]], [scalars_ts[1] * urms[1, 0], scalars_ts[-1] * urms[1, 0]], '--', c='k')
        # plt.legend()
        plt.title(r'log-log')

        plt.subplot(2, 3, 3)
        plt.plot(scalars_ts, np.real(urms[:, 0] - urms[0, 0]) / scalars_ts, label='(u_rms - u_rms(0))/t')
        # plt.loglog(ts, np.real(vrms[:, 0]), label='v_rms')
        # plt.loglog(ts, np.real(wrms[:, 0]), label='w_rms')
        # plt.xlim(xmax=20)
        # plt.loglog([ts[10], ts[-1]], [ts[10]*urms[10, 0], ts[-1]*urms[10, 0]], '--', c='k')
        plt.xlim((0, 2 * LOP_ts[LOP_t0_ind]))
        plt.ylim(ymin=0)
        plt.legend()
        plt.title(r'compensated')

        plt.subplot(2, 3, 5)
        plt.loglog(scalars_ts, np.real(bxrms[:, 0]), label='bx_rms')
        plt.loglog(scalars_ts, np.real(byrms[:, 0]), label='by_rms')
        plt.loglog(scalars_ts, np.real(bzrms[:, 0]), label='bz_rms')
        # plt.xlim(xmax=20)
        plt.loglog([scalars_ts[10], scalars_ts[-1]], [scalars_ts[10] * bxrms[10, 0], scalars_ts[-1] * bxrms[10, 0]], '--', c='k')
        # plt.loglog([scalars_ts[1], scalars_ts[-1]], [scalars_ts[1] * bxrms[1, 0], scalars_ts[-1] * bxrms[1, 0]], '--', c='k')
        # plt.legend(loc='upper right')
        # plt.legend()
        plt.xlabel(r'$t$')

        plt.subplot(2, 3, 6)
        plt.plot(scalars_ts, np.real(bxrms[:, 0] - bxrms[0, 0]) / scalars_ts, label='(bx_rms - bx_rms(0))/t')
        # plt.loglog(ts, np.real(byrms[:, 0]), label='by_rms')
        # plt.loglog(ts, np.real(bzrms[:, 0]), label='bz_rms')
        # plt.xlim(xmax=20)
        # plt.loglog([ts[10], ts[-1]], [ts[10]*bxrms[10, 0], ts[-1]*bxrms[10, 0]], '--', c='k')
        # plt.legend(loc='upper right')
        plt.legend()
        plt.xlim((0, 2 * LOP_ts[LOP_t0_ind]))
        plt.ylim(ymin=0)
        plt.xlabel(r'$t$')
        plt.tight_layout()
        pdf.savefig()
        plt.close()

    # stuff from plot_dEdt_terms_from_scan.py
    if True:
        plt.subplot(3, 1, 1)
        plt.axvline(LOP_ts[LOP_t0_ind], c='k')
        plt.plot(scalars_ts, KEs, label='KE')
        # plt.plot(scalars_ts[:-1], delta_KE, label='d(KE)/dt')
        plt.plot(scalars_ts, MEs, label='ME')
        plt.legend()
        plt.xlim((scalars_ts[0], scalars_ts[-1]))
        plt.ylim(ymin=0)

        plt.subplot(3, 1, 2)
        plt.axvline(LOP_ts[LOP_t0_ind], c='k')
        plt.axhline(0, c='k')
        # plt.plot(scalars_ts, KE, label='KE')
        plt.plot(scalars_ts[:-1], delta_KE, label='d(KE)/dt')
        plt.plot(scalars_ts, np.real(K_prod[:, 0]), label='KE shear prod')
        plt.plot(scalars_ts, np.real(EM_work[:, 0]), label='EM work')
        plt.legend(loc='upper right')
        plt.xlim((scalars_ts[0], scalars_ts[-1]))
        ymax = np.max(np.abs(np.append(delta_KE, [np.real(K_prod[:, 0]), np.real(EM_work[:, 0])])))
        plt.ylim((-ymax, ymax))

        plt.subplot(3, 1, 3)
        plt.axvline(LOP_ts[LOP_t0_ind], c='k')
        plt.axhline(0, c='k')
        # plt.plot(scalars_ts, ME, label='ME')
        plt.plot(scalars_ts[:-1], delta_ME, label='d(ME)/dt')
        plt.plot(scalars_ts, np.real(M_prod[:, 0]), label='ME shear prod')
        plt.plot(scalars_ts, np.real(-EM_work[:, 0]), label='-(EM work)')
        plt.legend(loc='upper right')
        plt.xlim((scalars_ts[0], scalars_ts[-1]))
        ymax = np.max(np.abs(np.append(delta_ME, [np.real(M_prod[:, 0]), np.real(-EM_work[:, 0])])))
        plt.ylim((-ymax, ymax))
        plt.tight_layout()
        pdf.savefig()
        plt.close()

    # Slices of the LOP and the evolved state, lifted from plot_snapshots_composite_from_scan.py
    if True:
        z_basis = de.Chebyshev('z', Nz, interval=(-0.5 * Lz, 0.5 * Lz))
        x_basis = de.Fourier('x', 4, interval=(0.0, 2.0 * np.pi / kx))
        y_basis = de.Fourier('y', 4, interval=(0.0, 2.0 * np.pi / ky))
        domain = de.Domain([x_basis, y_basis, z_basis], grid_dtype=np.float64)
        x1 = x_basis.grid(scale=1)
        y1 = y_basis.grid(scale=1)
        z1 = z_basis.grid(scale=1)
        x16 = x_basis.grid(scale=16)
        y16 = y_basis.grid(scale=16)
        z16 = z_basis.grid(scale=16)

        field = domain.new_field(name='field')
        scale = 1.0

        for i, field_data in enumerate([field_data_LOP, field_data_target]):
            print(i)
            fig = plt.figure(figsize=(scale * 12.0, 8.0))
            gs0 = gridspec.GridSpec(2, 3, figure=fig)
            gs0_subs = []
            # axes = []
            for fi, fname in enumerate(field_names):
                print(fi)
                gs0_subs += [gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[fi])]
                ax1 = fig.add_subplot(gs0_subs[-1][0, 0])

                field.set_scales(1, keep_data=False)
                field['c'][1, 1, :] = field_data[fname]
                field.set_scales(16, keep_data=True)
                field_out = field['g']

                pmesh1 = ax1.pcolormesh(x16, z16, field_out[:-1, 0, :-1].T, shading='flat', rasterized=True)
                fig.colorbar(mappable=pmesh1, ax=ax1)
                ax1.set_title(fname + ', $y=0$')
                if fi == 0 or fi == 3:
                    ax1.set_ylabel(r'$z$')
                ax1.set_xlabel(r'$x$')

            if i == 0:
                plt.suptitle(r'LOP')
            else:
                plt.suptitle(r'Evolved state, $t={:4.2f}$'.format(LOP_ts[LOP_t0_ind]))
            plt.tight_layout()
            pdf.savefig()
            plt.close()
