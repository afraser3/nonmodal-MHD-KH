# Reads the files made by scan_linear_optimals_3d.py and plots the perturbations
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import h5py
from dedalus import public as de
from pathlib import Path
from configparser import ConfigParser
# import glob

run_name = 'run_kx-log-scan-saveLOPs'
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

# with h5py.File('LOPs_3D/run_kx-log-scan-saveLOPs/kxi_0001/kyi_0000/gains_mats.h5', mode='r') as file:
# for gains_file in gains_files:
with PdfPages(str(basedir / 'evals_kxscan.pdf')) as pdf:
    for kydir in kydirs[::-1]:
        print('starting: ', kydir)
        gains_file = str(kydir / 'gains_mats.h5')
        with h5py.File(gains_file, mode='r') as file:
            kx = np.array(file['scales/kx'])
            ky = np.array(file['scales/ky'])
            evals = np.array(file['setup_mats/eigenvalues'])

        plt.plot(np.real(evals), np.imag(evals), '.')
        plt.title(r'$k_x={:1.2e}$'.format(kx))
        plt.tight_layout()
        # plt.show()
        # plt.savefig(str(kydir / 'LOP_tind{}.png'.format(tind)))
        # plt.savefig(str(kydir / 'LOP_tind{}.pdf'.format(tind)))
        pdf.savefig()
        plt.close()
