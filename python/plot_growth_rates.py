import numpy as np
from matplotlib import pyplot as plt
import h5py
from configparser import ConfigParser
from pathlib import Path

# run = "run_kx-lin-ky-lin_MA1.1_Re250_kxmin0.1_T200_Nz512-k600-saveLOPsAndEvolved"
run = "run_kx-log-ky-lin_MA1.1_Re250_kxmin1e-3_T200_Nz512-k600-saveLOPsAndEvolved"
# The .cfg file for this run
config_file = Path("LOPs_3D/config_files/{}.cfg".format(run))
# data dir that was saved to
basedir = Path("LOPs_3D") / config_file.stem
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

gains_file = "LOPs_3D/{}/gains.h5".format(run)
with h5py.File(gains_file, mode='r') as file:
    kx_meshgrid = np.array(file['ks/kx'])
    ky_meshgrid = np.array(file['ks/ky'])
    kxs = kx_meshgrid[0]
    kys = ky_meshgrid[:, 0]

evalues = np.zeros((len(kxs), len(kys)), dtype=np.complex128)

for kxi, kx in enumerate(kxs):
    for kyi, ky in enumerate(kys):
        with h5py.File("LOPs_3D/{}/kxi_{:04d}/kyi_{:04d}/gains_mats.h5".format(run, kxi, kyi), mode='r') as file:
            eigenvalues = np.array(file['setup_mats/eigenvalues'])
            indx = np.argmax(np.real(eigenvalues))
            evalues[kxi, kyi] = eigenvalues[indx]

plt.pcolormesh(kxs, kys, np.real(evalues))
plt.colorbar()
plt.contour(kxs, kys, np.real(evalues), levels=[0])
plt.show()