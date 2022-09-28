import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import h5py

run = "run_A"
gains_file = "LOPs_3D/{}/gains.h5".format(run)
with h5py.File(gains_file, mode='r') as file:
    kxs = np.array(file['ks/kx'])
    kys = np.array(file['ks/ky'])
    Gs = np.array(file['gains/Gs'])

plt.contourf(kxs, kys, np.amax(Gs, axis=2), norm=colors.LogNorm())
plt.colorbar()
plt.xlabel(r'$k_x$')
plt.ylabel(r'$k_y$')
# plt.show()
plt.savefig('LOPs_3D/{}/gains_kscan.png'.format(run))
plt.savefig('LOPs_3D/{}/gains_kscan.pdf'.format(run))
