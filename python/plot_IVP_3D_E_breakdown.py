import numpy as np
from matplotlib import pyplot as plt
import h5py
from pathlib import Path

test_ind = 46

Nz = 256
Reynolds = 50.0
Pm = 1.0
mReynolds = Pm * Reynolds
Lz = 10.0*np.pi
kx = 0.025
ky = 0.35
MA = 0.75
MA2 = MA**2.0
k = 300
dense = False

if dense:
    h5dir = Path("LOPs_3D/individual/Re{}/MA{}/Pm{}/kx{}/ky{}".format(Reynolds, MA, Pm, kx, ky))
    h5_fname = str(h5dir / 'Nz{}_dense.h5'.format(Nz))
    out_dir = h5dir / 'Nz{}_dense_ti{}'.format(Nz, test_ind)
    plot_name = str(out_dir / 'energies')
else:
    h5dir = Path("LOPs_3D/individual/Re{}/MA{}/Pm{}/kx{}/ky{}".format(Reynolds, MA, Pm, kx, ky))
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

KE = np.real(urms[:, 0])**2.0 + np.real(vrms[:, 0])**2.0 + np.real(wrms[:, 0])**2.0
ME = (np.real(bxrms[:, 0])**2.0 + np.real(byrms[:, 0])**2.0 + np.real(bzrms[:, 0])**2.0)/MA2

plt.subplot(2, 2, 1)
plt.plot(ts, np.real(urms[:, 0]), label='u_rms')
plt.plot(ts, np.real(vrms[:, 0]), label='v_rms')
plt.plot(ts, np.real(wrms[:, 0]), label='w_rms')
plt.legend()
plt.xlim((ts[0], ts[-1]))
plt.ylim(ymin=0)
plt.title('lin-lin')
plt.ylabel(r'$\mathbf{u}$ components')

plt.subplot(2, 2, 3)
plt.plot(ts, np.real(bxrms[:, 0]), label='bx_rms')
plt.plot(ts, np.real(byrms[:, 0]), label='by_rms')
plt.plot(ts, np.real(bzrms[:, 0]), label='bz_rms')
# plt.legend(loc='upper right')
plt.legend()
plt.xlim((ts[0], ts[-1]))
plt.ylim(ymin=0)
plt.ylabel(r'$\mathbf{b}$ components')
plt.xlabel(r'$t$')

plt.subplot(2, 2, 2)
plt.loglog(ts, np.real(urms[:, 0]), label='u_rms')
plt.loglog(ts, np.real(vrms[:, 0]), label='v_rms')
plt.loglog(ts, np.real(wrms[:, 0]), label='w_rms')
plt.legend()
plt.title(r'log-log')

plt.subplot(2, 2, 4)
plt.loglog(ts, np.real(bxrms[:, 0]), label='bx_rms')
plt.loglog(ts, np.real(byrms[:, 0]), label='by_rms')
plt.loglog(ts, np.real(bzrms[:, 0]), label='bz_rms')
# plt.legend(loc='upper right')
plt.legend()
plt.xlabel(r'$t$')

plt.savefig(str(out_dir) + '/energies_breakdown.png')
plt.savefig(str(out_dir) + '/energies_breakdown.pdf')
plt.show()
