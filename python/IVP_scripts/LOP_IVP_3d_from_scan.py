"""
Reads the files made by scan_linear_optimals_3d.py and does an IVP

Can either tell it a kx, ky, and t0, or just a t0 and it will find the largest-G kx/ky at that t0
"""

import numpy as np
from matplotlib import pyplot as plt
import h5py
from configparser import ConfigParser
from pathlib import Path
from dedalus import public as de
import logging
logger = logging.getLogger(__name__)

run = "run_kx-log-ky-lin_MA0.5_Re50_kxmin1e-3_T3tA_Nz256-k300-saveLOPs"
LOP_t0_ind = 11
kxi = None
kyi = None

if kxi is None and kyi is None:
    get_max_ks = True
else:
    get_max_ks = False

# The .cfg file for this run
config_file = Path("../LOPs_3D/config_files/{}.cfg".format(run))
# data dir that was saved to
basedir = Path("../LOPs_3D") / config_file.stem
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

t_start = params.getfloat('t_start')
t_stop = params.getfloat('t_stop')
t_num = params.getint('t_num')
t_endpoint = params.getboolean('t_endpoint')
times = np.linspace(t_start, t_stop, t_num, endpoint=t_endpoint)
save_LOPs_istep = params.getint('save_LOPs_istep')
LOP_ts = times[::save_LOPs_istep]

k = params.getint('k')
dense = params.getboolean('dense')
Nz = params.getint('Nz')
Lz = params.getfloat('Lz_factor')*np.pi

gains_file = "../LOPs_3D/{}/gains.h5".format(run)
with h5py.File(gains_file, mode='r') as file:
    kx_meshgrid = np.array(file['ks/kx'])
    ky_meshgrid = np.array(file['ks/ky'])
    kxs = kx_meshgrid[0]
    kys = ky_meshgrid[:, 0]
    Gs_global = np.array(file['gains/Gs'])
    Gs = Gs_global[:, :, ::save_LOPs_istep]
    # times = np.array(file['scales/ts'])

if get_max_ks:
    Gs_flat = Gs.reshape((-1, Gs.shape[-1]))  # flatten into a 2D array of shape (Gs.shape[0]*Gs.shape[1], Gs.shape[2])
    ind_flat = Gs_flat.argmax(axis=0)  # get the index where Gs_flat is maximized for each t
    kyind, kxind = np.unravel_index(ind_flat, Gs.shape[:2])  # unflatten/unravel those indices back into their kx, ky shape
    kxky_dirs = [basedir / 'kxi_{:04d}/kyi_{:04d}/'.format(kxind[i], kyind[i]) for i in range(len(Gs[0, 0]))]  # a list where each element points to the directory corresponding to the kx/ky that gives max G at each time
    kxky_dir = kxky_dirs[LOP_t0_ind]
    kxi = kxind[LOP_t0_ind]
    kyi = kyind[LOP_t0_ind]
    kx = kxs[kxi]
    ky = kys[kyi]
else:
    kxky_dir = basedir / 'kxi_{:04d}/kyi_{:04d}/'.format(kxi, kyi)
    kx = kxs[kxi]
    ky = kys[kyi]
if get_max_ks:
    out_dir = basedir / 'IVP_kxi{:04d}_kyi{:04d}_t0i{:03d}'.format(kxi, kyi, LOP_t0_ind)
else:
    out_dir = basedir / 'IVP_specific_kxi{:04d}_kyi{:04d}_t0i{:03d}'.format(kxi, kyi, LOP_t0_ind)
if not out_dir.exists():
    out_dir.mkdir(parents=False)

with h5py.File(str(kxky_dir / 'gains_mats.h5'), mode='r') as file:
    Gs_local = np.array(file['gains/Gs'])
    LOP = np.array(file['gains/LOPs'])[:, LOP_t0_ind]
    eigenvalues = np.array(file['setup_mats/eigenvalues'])
    indx = np.argmax(np.real(eigenvalues))

z_basis = de.Chebyshev('z', Nz, interval=(-0.5 * Lz, 0.5 * Lz))
domain = de.Domain([z_basis], grid_dtype=np.complex128)
probvars = ['u', 'v', 'w', 'p', 'Ax', 'Ay', 'Az', 'Phi', 'omega_x', 'omega_y', 'Bx', 'By']

problem = de.IVP(domain, variables=probvars, time='t')
z = domain.grid(0)
a = 1.0

nccU = domain.new_field(name='U')
nccU['g'] = np.tanh(z / a)
problem.parameters['U'] = nccU

nccDU = domain.new_field(name='DU')
nccDU['g'] = 1 / (a * np.cosh(z / a) ** 2)
problem.parameters['DU'] = nccDU

problem.parameters['MA2'] = MA2  # Alfven Mach number squared
problem.parameters['kx'] = kx
problem.parameters['ky'] = ky
problem.parameters['Re'] = Reynolds
problem.parameters['Rm'] = mReynolds
problem.substitutions['dx(A)'] = "1.0j*kx*A"
problem.substitutions['dy(A)'] = "1.0j*ky*A"

problem.substitutions['Bz'] = "dx(Ay) - dy(Ax)"
problem.substitutions['omega_z'] = "dx(v) - dy(u)"

# Coulumb gauge
problem.substitutions['Jx'] = "dy(Bz) - dz(By)"
problem.substitutions['Jy'] = "dz(Bx) - dx(Bz)"
problem.substitutions['Jz'] = "dx(By) - dy(Bx)"

problem.add_equation("dt(u) + U*dx(u) + w*DU + dx(p)          + (dy(omega_z) - dz(omega_y))/Re = 0")
problem.add_equation("dt(v) + U*dx(v)        + dy(p) - Jz/MA2 + (dz(omega_x) - dx(omega_z))/Re = 0")
problem.add_equation("dt(w) + U*dx(w)        + dz(p) + Jy/MA2 + (dx(omega_y) - dy(omega_x))/Re = 0")
problem.add_equation("dx(u) + dy(v) + dz(w) = 0")
problem.add_equation("dy(w) - dz(v) - omega_x = 0")
problem.add_equation("dz(u) - dx(w) - omega_y = 0")

problem.add_equation("dt(Ax) + dx(Phi)             + Jx/Rm = 0")
problem.add_equation("dt(Ay) + dy(Phi) - w + U*Bz  + Jy/Rm = 0")
problem.add_equation("dt(Az) + dz(Phi) + v - U*By  + Jz/Rm = 0")
problem.add_equation("dx(Ax) + dy(Ay) + dz(Az) = 0")
problem.add_equation("dy(Az) - dz(Ay) - Bx  = 0")
problem.add_equation("dz(Ax) - dx(Az) - By  = 0")

# no slip
problem.add_bc("left(u) = 0")
problem.add_bc("left(v) = 0")
problem.add_bc("left(w) = 0")
problem.add_bc("right(u) = 0")
problem.add_bc("right(v) = 0")
problem.add_bc("right(w) = 0")

# perfectly conducting boundaries
problem.add_bc("right(Phi) = 0")
problem.add_bc("right(Ax) = 0")
problem.add_bc("right(Ay) = 0")
problem.add_bc("left(Phi) = 0")
problem.add_bc("left(Ax) = 0")
problem.add_bc("left(Ay) = 0")

logger.info("done with BCs")

solver = problem.build_solver(de.timesteppers.RK443)
solver.stop_sim_time = t_stop
for i, var in enumerate(probvars):
    field = solver.state[var]
    field['c'] = LOP[int(i*len(LOP)/len(probvars)):int((i+1)*len(LOP)/len(probvars))]


def energy_norm(X1):
    u1 = X1['u']
    v1 = X1['v']
    w1 = X1['w']
    bx1 = X1['Bx']
    by1 = X1['By']
    bz1 = 1j * kx * X1['Ay'] - 1j * ky * X1['Ax']
    field = (np.abs(u1) ** 2.0 + np.abs(v1) ** 2.0 + np.abs(w1) ** 2.0 + (1.0 / MA2) * (np.abs(bx1) ** 2.0 + np.abs(by1) ** 2.0 + np.abs(bz1) ** 2.0)).evaluate().integrate()
    return field['g'][0]


def kinetic_energy(X1):
    u1 = X1['u']
    v1 = X1['v']
    w1 = X1['w']
    field = (np.abs(u1) ** 2.0 + np.abs(v1) ** 2.0 + np.abs(w1) ** 2.0).evaluate().integrate()
    return field['g'][0]


def magnetic_energy(X1):
    bx1 = X1['Bx']
    by1 = X1['By']
    bz1 = 1j * kx * X1['Ay'] - 1j * ky * X1['Ax']
    field = ((1.0 / MA2) * (np.abs(bx1) ** 2.0 + np.abs(by1) ** 2.0 + np.abs(bz1) ** 2.0)).evaluate().integrate()
    return field['g'][0]


scalars = solver.evaluator.add_file_handler(str(out_dir / 'scalars'), iter=1)
scalars.add_task('sqrt(integ(abs(u)**2))', layout='g', name='u_rms')
scalars.add_task('sqrt(integ(abs(v)**2))', layout='g', name='v_rms')
scalars.add_task('sqrt(integ(abs(w)**2))', layout='g', name='w_rms')
scalars.add_task('sqrt(integ(abs(Bx)**2))', layout='g', name='Bx_rms')
scalars.add_task('sqrt(integ(abs(By)**2))', layout='g', name='By_rms')
scalars.add_task('sqrt(integ(abs(Bz)**2))', layout='g', name='Bz_rms')
# scalars.add_task('integ(-w*u*DU)', layout='g', name='K_shear_prod')
# scalars.add_task('integ(Bz*Bx*DU)', layout='g', name='M_shear_prod')
# scalars.add_task('integ(dx(u*Bx + v*By + w*Bz))', layout='g', name='EM_work')
scalars.add_task('integ(-(w*conj(u) + conj(w)*u)*DU)', layout='g', name='K_shear_prod')
scalars.add_task('integ((Bz*conj(Bx) + conj(Bz)*Bx)*DU)/MA2', layout='g', name='M_shear_prod')
# scalars.add_task('integ(dx(u*conj(Bx) + conj(u)*Bx + v*conj(By) + conj(v)*By + w*conj(Bz) + conj(w)*Bz))', layout='g', name='EM_work')
scalars.add_task('-integ(Bx*conj(dx(u)) + conj(Bx)*dx(u) + By*conj(dx(v)) + conj(By)*dx(v) + Bz*conj(dx(w)) + conj(Bz)*dx(w))/MA2', layout='g', name='EM_work')

profiles = solver.evaluator.add_file_handler(str(out_dir / 'profiles'), iter=1)
profiles.add_task('u', layout='g', name='u')
profiles.add_task('v', layout='g', name='v')
profiles.add_task('w', layout='g', name='w')
profiles.add_task('omega_x', layout='g', name='omega_x')
profiles.add_task('omega_y', layout='g', name='omega_y')
profiles.add_task('omega_z', layout='g', name='omega_z')
profiles.add_task('Bx', layout='g', name='Bx')
profiles.add_task('By', layout='g', name='By')
profiles.add_task('Bz', layout='g', name='Bz')
profiles.add_task('Jx', layout='g', name='Jx')
profiles.add_task('Jy', layout='g', name='Jy')
profiles.add_task('Jz', layout='g', name='Jz')

profiles_c = solver.evaluator.add_file_handler(str(out_dir / 'profiles_c'), iter=1)
profiles_c.add_task('u', layout='c', name='u')
profiles_c.add_task('v', layout='c', name='v')
profiles_c.add_task('w', layout='c', name='w')
profiles_c.add_task('omega_x', layout='c', name='omega_x')
profiles_c.add_task('omega_y', layout='c', name='omega_y')
profiles_c.add_task('omega_z', layout='c', name='omega_z')
profiles_c.add_task('Bx', layout='c', name='Bx')
profiles_c.add_task('By', layout='c', name='By')
profiles_c.add_task('Bz', layout='c', name='Bz')
profiles_c.add_task('Jx', layout='c', name='Jx')
profiles_c.add_task('Jy', layout='c', name='Jy')
profiles_c.add_task('Jz', layout='c', name='Jz')

# should probably make dt min(1, 1/omega_Alfven) or something
# I think omega_Alfven is k_x/MA
# dt = 1e-1  # 2.0  # 2e-2
dt = min(1.0, MA/kx)
print('dt: ', dt)
ts = np.arange(0.0, t_stop + dt, dt)
energies = np.zeros_like(ts)
KEs = np.zeros_like(ts)
MEs = np.zeros_like(ts)
ti = 0
while solver.ok:
    energies[ti] = energy_norm(solver.state)
    KEs[ti] = kinetic_energy(solver.state)
    MEs[ti] = magnetic_energy(solver.state)
    solver.step(dt)
    if solver.iteration % 100 == 0:
        log_string = 'Iteration: {:5d}, '.format(solver.iteration)
        log_string += 'Time: {:8.3e}, '.format(solver.sim_time)
        log_string += 'Energy: {:8.3e}, '.format(energies[ti])
        logger.info(log_string)
    ti += 1
energies[-1] = energy_norm(solver.state)
KEs[-1] = kinetic_energy(solver.state)
MEs[-1] = magnetic_energy(solver.state)

plt.plot(times, np.amax(Gs_global, axis=(0, 1)), c='k', label=r'$G(t_0)$ maxed over kx/ky')
plt.plot(times, Gs_local, c='k', ls='--', label=r'$G(t_0)$ for this kx/ky')
plt.plot(LOP_ts[LOP_t0_ind], Gs_local[::save_LOPs_istep][LOP_t0_ind], 'x', c='k')
# plt.plot(ts, energies, label=r'$E(t, t_0={:4.2f})$'.format(ts_LOP[test_ind]))
# plt.plot(ts, KEs, label=r'$KE(t, t_0={:4.2f})$'.format(ts_LOP[test_ind]))
# plt.plot(ts, MEs, label=r'$ME(t, t_0={:4.2f})$'.format(ts_LOP[test_ind]))
plt.plot(ts, energies, label=r'$KE+ME$')
plt.plot(ts, KEs, label=r'$KE$')
plt.plot(ts, MEs, label=r'$ME$')
plt.xlabel(r'$t, t_0$')
plt.ylabel(r'$G, E$')
plt.xlim((ts[0], ts[-1]))
plt.legend()
plot_name = str(out_dir / 'energies')
plt.savefig(plot_name+'.png')
plt.savefig(plot_name+'.pdf')
plt.yscale('log')
plt.semilogy(ts, Gs_local[0]*np.exp(ts*2.0*np.max(np.real(eigenvalues[indx])))*Gs_local[-1]/np.exp(ts[-1]*2.0*np.max(np.real(eigenvalues[indx]))), '--', c='grey', label=r'eigenvalue')
plt.legend()
plt.ylim(ymin=0)
plt.tight_layout()
plt.savefig(plot_name+'_log.png')
plt.savefig(plot_name+'_log.pdf')
plt.show()
