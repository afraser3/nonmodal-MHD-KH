import numpy as np
from matplotlib import pyplot as plt
plt.style.use('prl.mplstyle')
# import matplotlib.colors as colors
# from matplotlib import gridspec as gridspec
import h5py
from configparser import ConfigParser
from pathlib import Path
from dedalus import public as de

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
run_dict_pleiades = {('big', 0.1): "run_kx-log-ky-lin_MA0.1_Re250_kxmin1e-3_T3tA_Nz512-k600-saveLOPsAndEvolved",
            ('big', 0.25): "run_kx-log-ky-lin_MA0.25_Re250_kxmin1e-3_T3tA_Nz512-k600-saveLOPsAndEvolved",
            ('big', 0.5): "run_kx-log-ky-lin_MA0.5_Re250_kxmin1e-3_T3tA_Nz512-k600-saveLOPsAndEvolved",
            ('big', 0.9): "run_kx-log-ky-lin_MA0.9_Re250_kxmin1e-3_T3tA_Nz512-k600-saveLOPsAndEvolved",
            ('big', 1): "run_kx-log-ky-lin_MA1_Re250_kxmin1e-3_T3tA_Nz512-k600-saveLOPsAndEvolved",
            ('big', 1.1, 5): "run_kx-log-ky-lin_MA1.1_Re250_kxmin1e-3_T5_Nz512-k600-saveLOPsAndEvolved",
            ('big', 1.1, 200): "run_kx-log-ky-lin_MA1.1_Re250_kxmin1e-3_T200_Nz512-k600-saveLOPsAndEvolved",
            ('big', 40): "run_kx-log-ky-lin_MA40_Re250_kxmin1e-3_T30_Nz512-k600-saveLOPsAndEvolved",
            ('small', 0.1): "run_kx-lin-ky-lin_MA0.1_Re250_kxmin0.1_T100_Nz512-k600-saveLOPsAndEvolved",
            ('small', 0.25): "run_kx-lin-ky-lin_MA0.25_Re250_kxmin0.1_T100_Nz512-k600-saveLOPsAndEvolved",
            ('small', 0.5): "run_kx-lin-ky-lin_MA0.5_Re250_kxmin0.1_kymax5_T100_Nz512-k600-saveLOPsAndEvolved",
            ('small', 1): "run_kx-lin-ky-lin_MA1_Re250_kxmin0.1_kymax5_Nz512-k600-saveLOPsAndEvolved",
            ('small', 1.1): "run_kx-lin-ky-lin_MA1.1_Re250_kxmin0.1_kymax5_T200_Nz512-k600-saveLOPsAndEvolved",
            ('small', 40): "run_kx-lin-ky-lin_MA40_Re250_kxmin0.1_kymax5_T100_Nz512-k600-saveLOPsAndEvolved"}


def get_Gs_ks_from_run(run):
    run_data = {}  # take everything we think we might want and throw it into a dict

    # Parse the .cfg filename passed to script
    config_file = Path("../LOPs_3D/config_files/{}.cfg".format(run))

    # Parse settings from .cfg file and save to variables
    runconfig = ConfigParser()
    runconfig.read(str(config_file))

    # data dir that was saved to
    basedir = Path("../LOPs_3D") / config_file.stem
    run_data['basedir'] = basedir

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
    save_LOPs_istep = params.getint('save_LOPs_istep')
    LOP_ts = times[::save_LOPs_istep]

    run_data['times'] = times
    run_data['LOP_ts'] = LOP_ts
    run_data['save_LOPs_istep'] = save_LOPs_istep

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
    run_data['Nz'] = Nz
    run_data['Lz'] = Lz
    # kxdirs = [basedir / 'kxi_{:04d}/'.format(kx_i) for kx_i, kx in enumerate(kxs)]

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
    run_data['kxky_dirs'] = [basedir / 'kxi_{:04d}/kyi_{:04d}/'.format(kxind[i], kyind[i]) for i in range(len(Gs[0, 0]))]  # a list where each element points to the directory corresponding to the kx/ky that gives max G at each time
    run_data['kxinds'] = kxind
    run_data['kyinds'] = kyind
    run_data['kx_opt'] = kxs[kxind]
    run_data['ky_opt'] = kys[kyind]
    return run_data


def get_LOP_from_run(run, t0ind):
    # need to use the previous function to figure out which kx/ky subdirectory to open for the given t0ind
    # then open it up, figure out what field we actually want to plot, use dedalus to convert it from 1D to 2D, and return it
    run_data = get_Gs_ks_from_run(run)
    gains_file = str(run_data['kxky_dirs'][::run_data['save_LOPs_istep']][t0ind] / 'gains_mats.h5')
    with h5py.File(gains_file, mode='r') as file:
        kx = np.array(file['scales/kx'])
        ky = np.array(file['scales/ky'])
        LOP = np.array(file['gains/LOPs'])[:, t0ind]
    return LOP, kx[()], ky[()], run_data


def LOP_to_3D_array(domain, LOP, field_name, kx, ky, scale=1):
    probvars = ['u', 'v', 'w', 'p', 'Ax', 'Ay', 'Az', 'Phi', 'omega_x', 'omega_y', 'Bx', 'By']
    field = domain.new_field(name=field_name)
    if field_name != 'Bz':
        field_ind = probvars.index(field_name)
        slice_start = int(field_ind * len(LOP) / len(probvars))
        slice_stop = int((field_ind + 1) * len(LOP) / len(probvars))
        if ky > 0:
            field['c'][1, 1, :] = LOP[slice_start:slice_stop]
        else:
            field['c'][1, 0, :] = LOP[slice_start:slice_stop]
    else:
        Ax = domain.new_field(name='Ax')
        Ay = domain.new_field(name='Ay')

        Ax_ind = probvars.index('Ax')
        slice_start = int(Ax_ind * len(LOP) / len(probvars))
        slice_stop = int((Ax_ind + 1) * len(LOP) / len(probvars))
        if ky > 0:
            Ax['c'][1, 1, :] = LOP[slice_start:slice_stop]
        else:
            Ax['c'][1, 0, :] = LOP[slice_start:slice_stop]

        Ay_ind = probvars.index('Ay')
        slice_start = int(Ay_ind * len(LOP) / len(probvars))
        slice_stop = int((Ay_ind + 1) * len(LOP) / len(probvars))
        if ky > 0:
            Ay['c'][1, 1, :] = LOP[slice_start:slice_stop]
        else:
            Ay['c'][1, 0, :] = LOP[slice_start:slice_stop]
        field['c'] = 1j * kx * Ay['c'] - 1j * ky * Ax['c']  # Bz = dx(Ay) - dy(Ax)
    field.set_scales(scale, keep_data=True)
    return field['g']


def get_IVP_scalar_data(run, LOP_t0ind):
    run_data = get_Gs_ks_from_run(run)
    IVP_data = {}
    t0ind = LOP_t0ind*run_data['save_LOPs_istep']
    # gains_file = str(run_data['kxky_dirs'][::run_data['save_LOPs_istep']][t0ind] / 'gains_mats.h5')
    kxi = run_data['kxinds'][t0ind]
    kyi = run_data['kyinds'][t0ind]
    IVP_dir = run_data['basedir'] / 'IVP_kxi{:04d}_kyi{:04d}_t0i{:03d}'.format(kxi, kyi, LOP_t0ind)
    with h5py.File(str(IVP_dir / 'scalars/scalars_s1/scalars_s1_p0.h5'), mode='r') as file:
        IVP_data['urms'] = np.array(file['tasks/u_rms'])
        IVP_data['vrms'] = np.array(file['tasks/v_rms'])
        IVP_data['wrms'] = np.array(file['tasks/w_rms'])
        IVP_data['scalars_ts'] = np.array(file['scales/sim_time'])
        IVP_data['bxrms'] = np.array(file['tasks/Bx_rms'])
        IVP_data['byrms'] = np.array(file['tasks/By_rms'])
        IVP_data['bzrms'] = np.array(file['tasks/Bz_rms'])
        IVP_data['K_prod'] = np.array(file['tasks/K_shear_prod'])
        IVP_data['M_prod'] = np.array(file['tasks/M_shear_prod'])
        IVP_data['EM_work'] = np.array(file['tasks/EM_work'])
    return IVP_data
