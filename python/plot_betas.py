# Reads the files made by scan_linear_optimals_3d.py and plots the perturbations
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import h5py
from dedalus import public as de
from pathlib import Path
from configparser import ConfigParser
# import glob

run_name = 'run_kx-log-scan_Nz256-dense-newscipy-saveLOPs-saveEvecs-saveBetas'
# run_name = 'test_Nz128-dense-saveLOPs-saveEvecs-saveBetas'
basedir = Path("LOPs_3D") / run_name
config_file = Path("LOPs_3D/config_files/"+run_name+'.cfg')
# gains_files = glob.glob(str(basedir)+'/kxi*/kyi*/gains_mats.h5')
# kydirs = glob.glob(str(basedir)+'/kxi*/kyi*/')

# Parse settings from .cfg file and save to variables
runconfig = ConfigParser()
runconfig.read(str(config_file))

# Parameters
params = runconfig['parameters']
MA = params.getfloat('MA')
MA2 = MA**2.0
Reynolds = params.getfloat('Reynolds')
Pm = params.getfloat('Pm')
mReynolds = Pm*Reynolds
kx_log = params.getboolean('kx_log')
kx_start = params.getfloat('kx_start')
kx_stop = params.getfloat('kx_stop')
kx_num = params.getint('kx_num')
kx_endpoint = params.getboolean('kx_endpoint')
ky_start = params.getfloat('ky_start')
ky_stop = params.getfloat('ky_stop')
ky_num = params.getint('ky_num')
ky_endpoint = params.getboolean('ky_endpoint')
kys = np.linspace(ky_start, ky_stop, ky_num, endpoint=ky_endpoint)  # the full set of kys to scan over
if kx_log:
    kxs = np.geomspace(kx_start, kx_stop, kx_num, endpoint=kx_endpoint)  # the full set of kxs to scan over
else:
    kxs = np.linspace(kx_start, kx_stop, kx_num, endpoint=kx_endpoint)  # the full set of kxs to scan over
kxdirs = [basedir / 'kxi_{:04d}/'.format(kx_i) for kx_i, kx in enumerate(kxs)]
kydirs = [kxdir / 'kyi_0000/' for kxdir in kxdirs]

t_start = params.getfloat('t_start')
t_stop = params.getfloat('t_stop')
t_num = params.getint('t_num')
t_endpoint = params.getboolean('t_endpoint')
ts = np.linspace(t_start, t_stop, t_num, endpoint=t_endpoint)

k = params.getint('k')
dense = params.getboolean('dense')
Nz = params.getint('Nz')
Lz = params.getfloat('Lz_factor')*np.pi
make_plots = params.getboolean('make_plots')
plot_modes = params.getboolean('plot_modes')
save_LOPs = params.getboolean('save_LOPs')
save_LOPs_istep = params.getint('save_LOPs_istep')


# with h5py.File(str(basedir / 'kxi_0000/kyi_0000/gains_mats.h5'), mode='r') as file:
# with h5py.File(str(basedir / 'kxi_0021/kyi_0000/gains_mats.h5'), mode='r') as file:
# for kydir in kydirs[::-1]:
if True:
    kydir = kydirs[-1]
    gains_file = str(kydir / 'gains_mats.h5')
    with h5py.File(gains_file, mode='r') as file:
        evals_high = np.array(file['setup_mats/eigenvalues_high'])
        evals = np.array(file['setup_mats/eigenvalues'])
        evals_low = np.array(file['setup_mats/eigenvalues_low'])
        # low_order = np.argsort(np.imag(evals_low))
        # low_order = np.argsort(np.abs(evals_low))

        left_evals = np.array(file['setup_mats/left_eigenvalues'])
        # left_order = np.argsort(-np.imag(left_evals))
        # left_order = np.argsort(np.abs(left_evals))
        # left_order = np.zeros_like(low_order)
        # for i in range(len(left_order)):
            # left_order[i] = np.argmin(np.abs(np.conj(left_evals) - evals_low[low_order[i]]))

        LOPs = np.array(file['gains/LOPs'])
        betas = np.array(file['gains/betas'])
        # betas_reorder = betas[left_order, :]
        evecs_low = np.array(file['setup_mats/eigenvectors_low'])
        # evecs_low_reorder = evecs_low[:, low_order]
        evecs = np.array(file['setup_mats/eigenvectors'])
        # modified_left_evecs = np.array(file['gains/modified_left_eigenvectors'])
        # modified_left_evecs_reorder = modified_left_evecs[:, left_order]
        # Need to normalize the modified left eigenvectors such that they are biorthoNORMAL, not just biorthogonal, to the right eigenvectors
        # modified_left_evecs_reorder_norm = modified_left_evecs_reorder / np.conjugate(np.diag(modified_left_evecs_reorder.T.conj() @ evecs_low_reorder))
        # betas = modified_left_evecs_reorder_norm.T.conj() @ LOPs
        kx = np.array(file['scales/kx'])
        z = np.array(file['scales/z'])
        # print(np.allclose(left_evals, np.conj(evals_low)))
        # print(kx, np.sum(np.isclose(left_evals[left_order], np.conj(evals_low[low_order]))))
        LOP = LOPs[:, -1]
        # reconstruction = np.matmul(evecs_low_reorder, betas[:, -1])
        reconstruction = np.matmul(evecs_low, betas[:, -1])
        # u_LOP = LOP[:int(len(LOP)/12)]
        # u_rec = reconstruction[:int(len(LOP)/12)]
        # print(kx, np.allclose(LOP, reconstruction), np.allclose(left_evals[left_order], np.conj(evals_low[low_order])))
        print(kx, np.allclose(LOP, reconstruction), np.allclose(left_evals, np.conj(evals_low)))


# LOP = LOPs[:, -1]
# print(np.shape(LOP))
# print(np.shape(betas))
# print(np.shape(evecs_low))
# reconstruction = np.matmul(evecs_low, betas[:, -1])
# u_LOP = LOP[:int(len(LOP)/12)]
# u_rec = reconstruction[:int(len(LOP)/12)]
# u_rec = u_rec*u_LOP[1]/u_rec[1]
# plt.plot(z, np.real(u_LOP), '-', c='C0')
# plt.plot(z, np.imag(u_LOP), '-', c='C1')
# plt.plot(z, np.real(u_rec), '--', c='C2')
# plt.plot(z, np.imag(u_rec), '--', c='C3')
# plt.xlim((-16, -15))
plt.plot(np.real(evals_low), np.imag(evals_low), '.')
plt.plot(np.real(evals), np.imag(evals), '+')
# plt.plot(np.real(evals_high), np.imag(evals_high), 'x')
plt.plot(np.real(left_evals), -np.imag(left_evals), 'x')
plt.show()


# for gains_file in gains_files:
# with PdfPages(str(basedir / 'evals_kxscan.pdf')) as pdf:
#     for kydir in kydirs[::-1]:
#         print('starting: ', kydir)
#         gains_file = str(kydir / 'gains_mats.h5')
#         with h5py.File(gains_file, mode='r') as file:
#             kx = np.array(file['scales/kx'])
#             ky = np.array(file['scales/ky'])
#             evals = np.array(file['setup_mats/eigenvalues'])
#
#         plt.plot(np.real(evals), np.imag(evals), '.')
#         plt.title(r'$k_x={:1.2e}$'.format(kx))
#         plt.tight_layout()
#         # plt.show()
#         # plt.savefig(str(kydir / 'LOP_tind{}.png'.format(tind)))
#         # plt.savefig(str(kydir / 'LOP_tind{}.pdf'.format(tind)))
#         pdf.savefig()
#         plt.close()
