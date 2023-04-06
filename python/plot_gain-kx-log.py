import h5py
from matplotlib import pyplot as plt
from pathlib import Path
import numpy as np
plt.style.use('prl')

basedir = Path("LOPs_3D/run_kx-log-scan-saveLOPs")

with h5py.File(basedir/'gains.h5', mode='r') as file:
    # print(np.array(file['EVs']))
    # print(np.array(file['ks']))
    # print(np.array(file['gains']))
    # print(np.array(file['scales']))
    kxs = np.array(file['ks/kx'])
    ts = np.array(file['scales/ts'])
    Gs = np.array(file['gains/Gs'])
    nEVs = np.array(file['EVs/nEVs'])

plt.subplot(2, 1, 1)
plt.semilogx(kxs[0], nEVs[0])
plt.ylabel('n good EVs')

plt.subplot(2, 1, 2)
plt.semilogx(kxs[0], Gs[0, :, 100])
plt.ylabel(r'$G(kx, t=50)$')

# plt.show()
plt.savefig(basedir/'nEVs_Gs-vs-kx.pdf')
plt.savefig(basedir/'nEVs_Gs-vs-kx.png')
