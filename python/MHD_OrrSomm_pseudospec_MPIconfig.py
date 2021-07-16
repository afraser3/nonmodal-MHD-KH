"""
Calculates the epsilon-pseudospectrum for a shear layer in 2D incompressible MHD.
Scans over MA then saves results to one hdf5 file per MA.

To run in parallel with 2 processes, e.g., do:
    $ mpiexec -n 2 python3 MHD_OrrSomm_pseudospec_MPIconfig.py config_files/test.cfg
The above line would save results to runs/test/ (creating directories as necessary)

Config file options:
    --MA_start          First Alfven Mach number (MA) in scan [default: 0.8]
    --MA_stop           Endpoint of interval for MA scan [default: 1.4]
    --MA_num            Number of points in MA scan [default: 4]
    --MA_endpoint       Whether to include MA_stop in the scan or not [default: True]

    --Reynolds          Reynolds number [default: 500.0]
    --Pm                Magnetic Prandtl number [default: 1.0]
    --kx                Horizontal wavenumber [default: 0.4]

    --k                 Dimensionality of subspace pseudospec is drawn from [default: 100]
    --Nz                Resolution in z [default: 512]
    --Lz_factor         Box size in z divided by pi [default: 10.0]

    --make_plots        Whether to plot pseudospectra [default: False]

"""

from eigentools import Eigenproblem
from dedalus import public as de
import numpy as np
import sys
import h5py
import logging
from configparser import ConfigParser
from pathlib import Path
from mpi4py import MPI
CW = MPI.COMM_WORLD
logger = logging.getLogger(__name__)

# Parse the .cfg filename passed to script
config_file = Path(sys.argv[-1])
logger.info("Running with config file {}".format(str(config_file)))

# Parse settings from .cfg file and save to variables
runconfig = ConfigParser()
runconfig.read(str(config_file))

# Parameters
params = runconfig['parameters']
MA_start = params.getfloat('MA_start')
MA_stop = params.getfloat('MA_stop')
MA_num = params.getint('MA_num')
MA_endpoint = params.getboolean('MA_endpoint')
Reynolds = params.getfloat('Reynolds')
Pm = params.getfloat('Pm')
mReynolds = Pm*Reynolds
kx = params.getfloat('kx')

k = params.getint('k')
Nz = params.getint('Nz')
Lz = params.getfloat('Lz_factor')*np.pi

make_plots = params.getboolean('make_plots')
if make_plots:
    from matplotlib import pyplot as plt

MA_global = np.linspace(MA_start, MA_stop, MA_num, endpoint=MA_endpoint)  # the full set of MAs to scan over
MA_local = MA_global[CW.rank::CW.size]  # the MAs that this processor will scan over
MA0 = MA_global[0]  # just a placeholder MA for setting Dedalus EVP up
if CW.rank == 0:
    logger.info('Scanning over these MAs: {}'.format(MA_global))
logger.debug('This MPI process will scan over these MAs: {}'.format(MA_local))

# data dir to save to: if "config_files/run1.cfg" is input, then "runs/run1/" becomes the output directory
datadir = Path("runs") / config_file.stem
plotdir = datadir / 'plots'
# filenames that this MPI process will save to
filenames_global = [datadir / 'pseudospec_{:03d}.h5'.format(ma_i) for ma_i, ma in enumerate(MA_global)]
filenames_local = filenames_global[CW.rank::CW.size]

if CW.rank == 0:  # only do this for the 0th MPI process
    if not datadir.exists():
        datadir.mkdir(parents=True)
    if make_plots:
        if not plotdir.exists():
            plotdir.mkdir(parents=True)

psize = 100  # num grid points (each axis) in complex frequency space over which pseudospectrum is calculated
freq_axis_bounds = [-1.0 * kx / 0.4, 1.0 * kx / 0.4]  # the reasonable/helpful bounds seem to scale roughly with kx
growth_axis_bounds = [-1.0 * kx / 0.4, 0.2]  # 0.2 is sometimes too small to see full extent of epsilon=0.1
real_points = np.linspace(freq_axis_bounds[0], freq_axis_bounds[1], psize)
imag_points = np.linspace(growth_axis_bounds[0], growth_axis_bounds[1], psize)


def energy_norm_general(X1, X2, MA):  # because the norm depends on MA, need to include MA as an argument
    u1 = X1['phi_z']
    w1 = -1.0j * kx * X1['phi']
    bx1 = X1['psi_z']
    bz1 = -1.0j * kx * X1['psi']
    u2 = X2['phi_z']
    w2 = -1.0j * kx * X2['phi']
    bx2 = X2['psi_z']
    bz2 = -1.0j * kx * X2['psi']

    field = (np.conj(u1)*u2 + np.conj(w1)*w2 + (1.0/MA**2.0)*(np.conj(bx1)*bx2 + np.conj(bz1)*bz2)).evaluate().integrate()
    return field['g'][0]


if Reynolds == np.inf:
    inviscid = True
else:
    inviscid = False
if mReynolds == np.inf:
    ideal = True
else:
    ideal = False
z_basis = de.Chebyshev('z', Nz, interval=(-0.5 * Lz, 0.5 * Lz))
domain = de.Domain([z_basis], grid_dtype=np.complex128, comm=MPI.COMM_SELF)
probvars = ['phi', 'phi_z', 'psi', 'psi_z']
if not inviscid:
    probvars.append('uz')
    probvars.append('uzz')
else:  # really haven't tested this! May have neglected to add the right substitutions/equations
    probvars.append('zeta')
problem = de.EVP(domain, variables=probvars, eigenvalue='omega')
z = domain.grid(0)
a = 1.0

nccU = domain.new_field(name='U')
nccU['g'] = np.tanh(z / a)
problem.parameters['U'] = nccU

nccDU = domain.new_field(name='DU')
nccDU['g'] = 1 / (a * np.cosh(z / a) ** 2)
problem.parameters['DU'] = nccDU

problem.parameters['MA2'] = MA0 ** 2.0  # Alfven Mach number squared
problem.parameters['kx'] = kx
problem.parameters['Lz'] = Lz
if not inviscid:
    problem.parameters['Re'] = Reynolds
if not ideal:
    problem.parameters['Rm'] = mReynolds
problem.substitutions['dx(A)'] = "1.0j*kx*A"
problem.substitutions['dt(A)'] = "-1.0j*omega*A"
problem.substitutions['u'] = 'phi_z'
problem.substitutions['w'] = "-dx(phi)"
problem.substitutions['Bz'] = "-dx(psi)"
problem.substitutions['Bx'] = "psi_z"
problem.substitutions['Jz'] = "dx(dx(psi)) + dz(Bx)"
if not inviscid:
    problem.substitutions['zeta'] = "dx(dx(phi)) + uz"
    problem.substitutions['zetaz'] = "dx(dx(u)) + uzz"
if not inviscid:
    problem.add_equation(
        "dt(zeta) - 1/Re*(dx(dx(zeta)) + dz(zetaz)) + U*dx(zeta) - dx(phi)*dz(DU) - 1/MA2*(dx(Jz)) = 0")
    problem.add_equation("uz - dz(u) = 0")
    problem.add_equation("uzz - dz(uz) = 0")
else:
    problem.add_equation("dt(zeta) + U*dx(zeta) - dx(phi)*dz(DU) - 1/MA2*(dx(Jz)) = 0")
    problem.add_equation("zeta - dx(dx(phi)) - dz(phi_z) = 0")
if not ideal:
    problem.add_equation("dt(psi) - 1/Rm*Jz + U*dx(psi) - dx(phi) = 0")
else:  # haven't tested this either
    problem.add_equation("dt(psi) + U*dx(psi) - dx(phi) = 0")
problem.add_equation("phi_z - dz(phi) = 0")
problem.add_equation("psi_z - dz(psi) = 0")

problem.add_bc("left(w) = 0")
problem.add_bc("right(w) = 0")
if not inviscid:
    problem.add_bc("left(dz(w)) = 0")  # no-slip, co-moving (with base flow U_0) boundaries
    problem.add_bc("right(dz(w)) = 0")  # (does free-slip make more sense?)
problem.add_bc("right(Bz) = 0")  # perfectly conducting boundaries
problem.add_bc("left(Bz) = 0")

logger.info("done with BCs")

EP = Eigenproblem(problem, grow_func=lambda x: x.imag, freq_func=lambda x: x.real)

for ma_ind, ma in enumerate(MA_local):
    logger.info('computing pseudospectrum for MA={}'.format(ma))
    EP._set_parameters({'MA2': ma ** 2.0})
    EP.calc_ps(k, (real_points, imag_points), inner_product=lambda x1, x2: energy_norm_general(x1, x2, ma))
    with h5py.File(str(filenames_local[ma_ind]), 'w-') as file:
        evalues_grp = file.create_group('evalues')
        evecs_grp = file.create_group('evecs')
        pseudospec_grp = file.create_group('pseudospec')

        evalues = evalues_grp.create_dataset('evalues', data=EP.evalues)
        evalues_low = evalues_grp.create_dataset('evalues_low', data=EP.evalues_low)
        evalues_high = evalues_grp.create_dataset('evalues_high', data=EP.evalues_high)

        evectors = evecs_grp.create_dataset('evectors', data=EP.solver.eigenvectors)

        ps_real = pseudospec_grp.create_dataset('ps_real', data=EP.ps_real)
        ps_imag = pseudospec_grp.create_dataset('ps_imag', data=EP.ps_imag)
        pseudospectrum = pseudospec_grp.create_dataset('pseudospectrum', data=EP.pseudospectrum)
    if make_plots:
        plt.plot(np.real(EP.evalues), np.imag(EP.evalues), '.', c='k')
        if np.min(EP.pseudospectrum) < 0.1:
            plt.contour(EP.ps_real, EP.ps_imag, np.log10(EP.pseudospectrum), levels=np.arange(-8, 0))
            plt.colorbar(label=r'$\log_{10} (\epsilon)$')
        plt.ylim((EP.ps_imag[0], EP.ps_imag[-1]))
        plt.xlim((EP.ps_real[0], EP.ps_real[-1]))
        plt.axhline(0, color='k', alpha=0.2)
        plt.xlabel('real (oscillating) frequency')
        plt.ylabel('growth rate; max = {:.4e}'.format(np.max(np.imag(EP.evalues))))
        plt.title(r'(kx, MA, Re, Pm, k) = ({}, {:4.4f}, {}, {}, {})'.format(kx, ma, Reynolds, Pm, k))
        for ext in ['.pdf', '.png']:
            plt.savefig(str(plotdir / str(filenames_local[ma_ind].stem)) + ext)
        plt.close()
