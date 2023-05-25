import numpy as np
from matplotlib import pyplot as plt
import h5py
from pathlib import Path

test_ind = 287

Nz = 512
Reynolds = 250.0
Pm = 1.0
mReynolds = Pm * Reynolds
Lz = 10.0*np.pi
kx = 0.025
ky = 0.375
MA = 0.25
MA2 = MA**2.0
k = 300
dense = False

if dense:
    h5dir = Path("../LOPs_3D/individual/Re{}/MA{}/Pm{}/kx{}/ky{}".format(Reynolds, MA, Pm, kx, ky))
    h5_fname = str(h5dir / 'Nz{}_dense.h5'.format(Nz))
    out_dir = h5dir / 'Nz{}_dense_ti{}'.format(Nz, test_ind)
    plot_name = str(out_dir / 'energies')
else:
    h5dir = Path("../LOPs_3D/individual/Re{}/MA{}/Pm{}/kx{}/ky{}".format(Reynolds, MA, Pm, kx, ky))
    h5_fname = str(h5dir / 'Nz{}_k{}.h5'.format(Nz, k))
    out_dir = h5dir / 'Nz{}_k{}_ti{}'.format(Nz, k, test_ind)
    plot_name = str(out_dir / 'energies')

with h5py.File(str(out_dir) + '/scalars/scalars_s1/scalars_s1_p0.h5', mode='r') as file:
    urms = np.array(file['tasks/u_rms'])
    vrms = np.array(file['tasks/v_rms'])
    wrms = np.array(file['tasks/w_rms'])
    ts = np.array(file['scales/sim_time'])
    bxrms = np.array(file['tasks/Bx_rms'])
    byrms = np.array(file['tasks/By_rms'])
    bzrms = np.array(file['tasks/Bz_rms'])
    K_prod = np.array(file['tasks/K_shear_prod'])
    M_prod = np.array(file['tasks/M_shear_prod'])
    EM_work = np.array(file['tasks/EM_work'])

KE = np.real(urms[:, 0])**2.0 + np.real(vrms[:, 0])**2.0 + np.real(wrms[:, 0])**2.0
ME = (np.real(bxrms[:, 0])**2.0 + np.real(byrms[:, 0])**2.0 + np.real(bzrms[:, 0])**2.0)/MA2
delta_KE = np.diff(KE)/ts[1]
delta_ME = np.diff(ME)/ts[1]
plt.subplot(3, 1, 1)
plt.plot(ts, KE, label='KE')
# plt.plot(ts[:-1], delta_KE, label='d(KE)/dt')
plt.plot(ts, ME, label='ME')
# plt.plot(ts[:-1], delta_ME, label='d(ME)/dt')
# plt.plot(ts, np.real(K_prod[:, 0]), label='KE shear prod')
# plt.plot(ts, np.real(M_prod[:, 0]), label='ME shear prod')
# plt.plot(ts, np.real(EM_work[:, 0]), label='EM work')
plt.legend()
plt.xlim((ts[0], ts[-1]))
plt.ylim(ymin=0)
plt.subplot(3, 1, 2)
# plt.plot(ts, KE, label='KE')
plt.plot(ts[:-1], delta_KE, label='d(KE)/dt')
plt.plot(ts, np.real(K_prod[:, 0]), label='KE shear prod')
plt.plot(ts, np.real(EM_work[:, 0]), label='EM work')
plt.legend(loc='upper right')
plt.xlim((ts[0], ts[-1]))
plt.subplot(3, 1, 3)
# plt.plot(ts, ME, label='ME')
plt.plot(ts[:-1], delta_ME, label='d(ME)/dt')
plt.plot(ts, np.real(M_prod[:, 0]), label='ME shear prod')
plt.plot(ts, np.real(-EM_work[:, 0]), label='-(EM work)')
plt.legend(loc='upper right')
plt.xlim((ts[0], ts[-1]))
plt.savefig(str(out_dir) + '/energies_prods.png')
plt.savefig(str(out_dir) + '/energies_prods.pdf')
plt.show()
