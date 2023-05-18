"""
Wrote this for plotting nEVs and Gs from scans I did where I had just one ky and a bunch of kx's spaced out evenly in
log space, so note this will break in various places if multiple ky's were run!
"""

import h5py
from matplotlib import pyplot as plt
from pathlib import Path
import numpy as np
import glob
from configparser import ConfigParser
import re
plt.style.use('prl')


def natural_sort(iterable, reverse=False):
    """
    Sort alphanumeric strings naturally, i.e. with "1" before "10".
    Based on http://stackoverflow.com/a/4836734.

    Actually, just copy+pasted from dedalus
    """

    convert = lambda sub: int(sub) if sub.isdigit() else sub.lower()
    key = lambda item: [convert(sub) for sub in re.split('([0-9]+)', str(item))]

    return sorted(iterable, key=key, reverse=reverse)


def get_data_from_bad_kxs(run_name, runconfig):
    # Parse settings from .cfg file and save to variables
    # If I wrote this better I wouldn't have to pass runconfig as well as run_name...
    # config_file = Path("LOPs_3D/config_files/" + run_name + '.cfg')
    # runconfig = ConfigParser()
    # runconfig.read(str(config_file))
    # Parameters
    params = runconfig['parameters']

    if params.getboolean('kx_log'):
        kx_global = np.geomspace(params.getfloat('kx_start'), params.getfloat('kx_stop'), params.getint('kx_num'), endpoint=params.getboolean('kx_endpoint'))  # the full set of kxs that was attempted
    else:
        kx_global = np.linspace(params.getfloat('kx_start'), params.getfloat('kx_stop'), params.getint('kx_num'), endpoint=params.getboolean('kx_endpoint'))
    if params.getboolean('add_kx_zero') and kx_global[0] != 0:
        kx_global = np.insert(kx_global, 0, 0.0)
    kxs = np.empty((1, len(kx_global)), dtype=np.float64)
    kxs[0] = kx_global
    gainsfiles = natural_sort(glob.glob(str(Path("LOPs_3D") / run_name) + '/kxi*/kyi*/gains_mats.h5'))
    # gainfiles is a list of strings, where each string is the path to one of the gains_mats.h5 files in this run. What
    # I want to do is go through each item in that list and pick out the "kxi" in the path to identify which kx's
    # successfully generated a gains_mats.h5 file. The below does that through a bunch of tedious string.split() stuff
    good_kxinds = [int(file.split('/')[2].split('_')[1]) for file in gainsfiles]
    for i, kxi in enumerate(good_kxinds):
        with h5py.File(gainsfiles[i], mode='r') as file:
            try:
                gains = np.array(file['gains/Gs'])
            except KeyError:
                print(kxi)
                continue
            num_evs = len(np.array(file['setup_mats/eigenvalues']))
            if i == 0:
                Gs = np.empty((1, len(kx_global), len(gains)), dtype=np.float64)
                Gs[:] = np.nan
                nEVs = np.empty((1, len(kx_global)), dtype=np.float64)
                nEVs[:] = np.nan
                ts = np.array(file['scales/ts'])
            Gs[0, kxi] = gains
            nEVs[0, kxi] = num_evs
    return kxs, ts, Gs, nEVs


# basedir = Path("LOPs_3D/run_kx-log-scan-MA40_Nz256-dense-newscipy-saveLOPs-saveEvecs-saveBetas")
# basedir2 = Path("LOPs_3D/run_kx-log-scan_Nz256-dense-newscipy-saveLOPs-saveEvecs-saveBetas")
# basedir = Path("LOPs_3D/run_kx-log-scan_Nz512-k600")
run_name = "run_kx-log-scan_Nz512-dense-Legendre-saveLOPs-saveEvecs-saveBetas"
compare = True
run_name2 = "run_kx-log-scan_Nz512-dense-newscipy-saveLOPs-saveEvecs-saveBetas"

basedir = Path("LOPs_3D") / run_name

if compare:
    basedir2 = Path("LOPs_3D") / run_name2

bad_kxs = True  # workaround in case the main gains.h5 file doesn't exist -- go through all the individual kxi folders

if bad_kxs:
    config_file = Path("LOPs_3D/config_files/" + run_name + '.cfg')
    runconfig = ConfigParser()
    runconfig.read(str(config_file))
    kxs, ts, Gs, nEVs = get_data_from_bad_kxs(run_name, runconfig)
    if compare:
        config_file2 = Path("LOPs_3D/config_files/" + run_name2 + '.cfg')
        runconfig2 = ConfigParser()
        runconfig2.read(str(config_file2))
        kxs2, ts2, Gs2, nEVs2 = get_data_from_bad_kxs(run_name2, runconfig2)


else:
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

plt.figure(figsize=(5, 10))
plt.subplot(4, 1, 1)
plt.semilogx(kxs[0], nEVs[0], c='k')
if compare:
    plt.semilogx(kxs2[0], nEVs2[0], c='grey', ls='--')
plt.title('n good EVs')

colors = ['C0', 'C1', 'C2']
colors2 = ['C1', 'C2', 'C0']
tinds = [10, 33, 100]
for i, color in enumerate(colors):
    plt.subplot(4, 1, 2+i)
    plt.semilogx(kxs[0], Gs[0, :, tinds[i]], c=color)
    if compare:
        plt.semilogx(kxs2[0], Gs2[0, :, tinds[i]], c=colors2[i], ls='--')
    plt.title(r'$G(kx, t={:.2e})$'.format(ts[tinds[i]]))
plt.xlabel(r'$k_x$')

plt.tight_layout()
# plt.show()
plt.savefig(basedir/'nEVs_Gs-vs-kx.pdf')
plt.savefig(basedir/'nEVs_Gs-vs-kx.png')
