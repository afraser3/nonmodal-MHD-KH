import h5py
from matplotlib import pyplot as plt
from pathlib import Path
import numpy as np
plt.style.use('prl')

basedir = Path("LOPs_3D/run_kx-log-scan-MA40_Nz256-dense-newscipy-saveLOPs-saveEvecs-saveBetas")
# basedir2 = Path("LOPs_3D/run_kx-log-scan_Nz256-dense-newscipy-saveLOPs-saveEvecs-saveBetas")
# basedir = Path("LOPs_3D/run_kx-log-scan_Nz512-k600")
compare = False
basedir2 = Path("LOPs_3D/run_kx-log-scan-saveLOPs")

with h5py.File(basedir/'gains.h5', mode='r') as file:
    kxs = np.array(file['ks/kx'])
    ts = np.array(file['scales/ts'])
    Gs = np.array(file['gains/Gs'])
    nEVs = np.array(file['EVs/nEVs'])

if compare:
    with h5py.File(basedir2/'gains.h5', mode='r') as file:
        kxs2 = np.array(file['ks/kx'])
        ts2 = np.array(file['scales/ts'])
        Gs2 = np.array(file['gains/Gs'])
        nEVs2 = np.array(file['EVs/nEVs'])

plt.figure(figsize=(5,10))
plt.subplot(4, 1, 1)
plt.semilogx(kxs[0], nEVs[0], c='k')
if compare:
    plt.semilogx(kxs2[0], nEVs2[0], c='k', ls='--')
plt.title('n good EVs')

colors = ['C0', 'C1', 'C2']
tinds = [10, 33, 100]
for i, color in enumerate(colors):
    plt.subplot(4, 1, 2+i)
    plt.semilogx(kxs[0], Gs[0, :, tinds[i]], c=color)
    if compare:
        plt.semilogx(kxs2[0], Gs2[0, :, tinds[i]], c=color, ls='--')
    plt.title(r'$G(kx, t={:.2e})$'.format(ts[tinds[i]]))
plt.xlabel(r'$k_x$')

plt.tight_layout()
# plt.show()
plt.savefig(basedir/'nEVs_Gs-vs-kx.pdf')
plt.savefig(basedir/'nEVs_Gs-vs-kx.png')
