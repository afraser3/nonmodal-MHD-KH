"""
Calculates the epsilon-pseudospectrum for a shear layer in 2D incompressible MHD.
Scans over MA then saves results to one hdf5 file per MA.

TODO: logger.info only shows CW.rank==0 stuff. What's the best way to output the MA_local of other processes too?
TODO: might as well make the plots with this script too, right? Maybe with a flag put in the .cfg file
TODO: how can I put logger settings in the .cfg file? Like whether to log warnings, info, debug, etc.

To run in parallel with 2 processes, e.g., do:
    $ mpiexec -n 2 python3 MHD_OrrSomm_pseudospec_MPIconfig.py config_files/test.cfg
The usage below is for running in serial (haven't merged with serial case because I'm new to docopt)

Usage:
    MHD_OrrSomm_pseudospec_MPIconfig.py [options]
    MHD_OrrSomm_pseudospec_MPIconfig.py <config> [options]


Options:
    --MA_start=<MA_start>       First Alfven Mach number (MA) in scan [default: 0.8]
    --MA_stop=<MA_stop>         Endpoint of interval for MA scan [default: 1.4]
    --MA_num=<MA_num>           Number of points in MA scan [default: 4]
    --MA_endpoint=<MA_endpoint> Whether to include MA_stop in the scan or not [default: True]

    --Reynolds=<Reynolds>       Reynolds number [default: 500.0]
    --Pm=<Pm>                   Magnetic Prandtl number [default: 1.0]
    --kx=<kx>                   Horizontal wavenumber [default: 0.4]

    --k=<k>                     Dimensionality of subspace pseudospec is drawn from [default: 100]
    --Nz=<Nz>                   Resolution in z [default: 512]
    --Lz_factor=<Lz>            Box size in z divided by pi [default: 10.0]

    --make_plots=<make_plots>   Whether to plot pseudospectra [default: False]

    --root_dir=<dir>            Root directory for output [default: ./]
    --label=<label>             Add label to directory name

"""

from eigentools import Eigenproblem
from dedalus import public as de
import numpy as np
import os
import h5py
import logging
from configparser import ConfigParser
from pathlib import Path
from docopt import docopt
from mpi4py import MPI
CW = MPI.COMM_WORLD
logger = logging.getLogger(__name__)

args = docopt(__doc__)
if args['<config>'] is not None:
    config_file = Path(args['<config>'])
    config = ConfigParser()
    config.read(str(config_file))
    for n, v in config.items('parameters'):
        for k in args.keys():
            if k.split('--')[-1].lower() == n:  # what does this do?
                if v == 'true':
                    v = True
                args[k] = v

MA_start = float(args['--MA_start'])
MA_stop = float(args['--MA_stop'])
MA_num = int(args['--MA_num'])
MA_endpoint = bool(args['--MA_endpoint'])
Reynolds = float(args['--Reynolds'])
Pm = float(args['--Pm'])
mReynolds = Pm*Reynolds
kx = float(args['--kx'])

k = int(args['--k'])
Nz = int(args['--Nz'])
Lz = float(args['--Lz_factor'])*np.pi

make_plots = bool(args['--make_plots'])
if make_plots:
    from matplotlib import pyplot as plt

MA_global = np.linspace(MA_start, MA_stop, MA_num, endpoint=MA_endpoint)  # the full set of MAs to scan over
MA_local = MA_global[CW.rank::CW.size]  # the MAs that this processor will scan over
if CW.rank == 0:
    logger.info('Scanning over these MAs: {}'.format(MA_global))
logger.debug('This MPI process will scan over these MAs: {}'.format(MA_local))

MA0 = MA_global[0]

psize = 100  # num grid points (each axis) in complex frequency space over which pseudospectrum is calculated
freq_axis_bounds = [-1.0 * kx / 0.4, 1.0 * kx / 0.4]  # the reasonable/helpful bounds seem to scale roughly with kx
growth_axis_bounds = [-1.0 * kx / 0.4, 0.2]  # 0.2 is sometimes too small to see full extent of epsilon=0.1
real_points = np.linspace(freq_axis_bounds[0], freq_axis_bounds[1], psize)
imag_points = np.linspace(growth_axis_bounds[0], growth_axis_bounds[1], psize)

base_path = 'saved_spectra/Pm{}/Re{}/'.format(Pm, Reynolds)
# filepaths = ['saved_spectra/Pm{}/Re{}/MA{}/kx{}/'.format(Pm, Reynolds, ma, kx) for ma in MA_global]
filepaths_local = ['saved_spectra/Pm{}/Re{}/MA{}/kx{}/'.format(Pm, Reynolds, ma, kx) for ma in MA_local]
filename = 'pseudospec_Nz{}_Lz{}pi_k{}.h5'.format(Nz, Lz/np.pi, k)

# ugly hack: skip_flag_local includes the list of MAs that should be skipped because those files already exist
# from a previous calculation (so if you scan over 200 MAs but 1 already has been done, no need to raise a whole error,
# just skip that one MA)
skip_flag_local = np.zeros_like(MA_local, dtype=bool)

if CW.rank == 0:  # only do this for the 0th MPI process
    try:
        os.makedirs(base_path)
    except FileExistsError:
        pass

for fi, filepath in enumerate(filepaths_local):
    try:
        os.makedirs(filepath)
    except FileExistsError:
        pass
    if filename in os.listdir(filepath):  # If a file already exists at this MA...
        skip_flag_local[fi] = True  # ...then note that for later on so we can skip it
        # raise FileExistsError  # previously I had been raising an error instead


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
# problem.substitutions['dt(A)'] = "1.0j*omega*A"
# problem.substitutions['dt(A)'] = "sigma*A"  # was having trouble with EP.calc_ps using omega
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

print("done with BCs")

EP = Eigenproblem(problem, grow_func=lambda x: x.imag, freq_func=lambda x: x.real)

for ma_ind, ma in enumerate(MA_local):
    if skip_flag_local[ma_ind]:
        logger.info('skipping calculation for MA={}'.format(ma))
    else:
        logger.info('computing pseudospectrum for MA={}'.format(ma))
        problem.namespace['MA2'].value = ma ** 2.0
        EP.calc_ps(k, (real_points, imag_points), inner_product=lambda x1, x2: energy_norm_general(x1, x2, ma))
        with h5py.File(filepaths_local[ma_ind]+filename, 'w-') as file:
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
            plotname = 'pseudospec_Nz{}_Lz{}pi_k{}.pdf'.format(Nz, Lz / np.pi, k)
            plt.plot(np.real(EP.evalues), np.imag(EP.evalues), '.', c='k')
            plt.contour(EP.ps_real, EP.ps_imag, np.log10(EP.pseudospectrum), levels=np.arange(-8, 0))
            plt.colorbar(label=r'$\log_{10} (\epsilon)$')
            plt.ylim((EP.ps_imag[0], EP.ps_imag[-1]))
            plt.xlim((EP.ps_real[0], EP.ps_real[-1]))
            plt.axhline(0, color='k', alpha=0.2)
            plt.xlabel('real (oscillating) frequency')
            plt.ylabel('growth rate')
            plt.title(r'(kx, MA, Re, Pm, k) = ({}, {}, {}, {}, {})'.format(kx, ma, Reynolds, Pm, k))
            plt.savefig(filepaths_local[ma_ind]+plotname)
            plt.close()
