# Reads the files made by LOP_save_mats.py and plots the perturbations
import numpy as np
from matplotlib import pyplot as plt
import h5py
from pathlib import Path
from dedalus import public as de
import logging
logger = logging.getLogger(__name__)

test_ind = 50

Nz = 256
Reynolds = 50.0
Pm = 1.0
mReynolds = Pm * Reynolds
Lz = 10.0*np.pi
kx = 0.025  # 0.2
ky = 0.35  # 0.025
MA = 1.0
MA2 = MA**2.0
k = 300
dense = False
reject_spurious = True

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
if not out_dir.exists():
    out_dir.mkdir(parents=False)

with h5py.File(h5_fname, mode='r') as file:
    eigenvalues = np.array(file['setup_mats/eigenvalues'])
    indx = np.argmax(np.real(eigenvalues))
    evec = np.array(file['setup_mats/eigenvectors'][:, indx])
    evec_phase = evec[0] / np.abs(evec[0])
    evec = evec / evec_phase
    F = np.array(file['setup_mats/F'])
    Finv = np.linalg.inv(F)
    # F_KE = np.array(file['setup_mats/F_KE'])
    # F_ME = np.array(file['setup_mats/F_ME'])
    R = np.array(file['setup_mats/R'])
    Q = np.array(file['setup_mats/Q'])
    pre_right = np.array(file['setup_mats/pre_right'])
    ts_LOP = np.array(file['scales/ts'])
    Nz = np.array(file['scales/Nz'])[()]
    Lz = np.array(file['scales/Lz'])[()]
    kx = np.array(file['scales/kx'])[()]
    # ky = np.array(file['scales/ky'])[()]
    LOPs = []
    # Gs = np.zeros_like(ts_LOP)
    Gs = np.array(file['energies/Gs'])
    # E_test = np.zeros_like(ts_LOP)
    # KE_test = np.zeros_like(ts_LOP)
    # ME_test = np.zeros_like(ts_LOP)
    # v1_test = np.array(file['ti_{}/Vh'.format(test_ind)]).conj()[0]
    v1_test = np.array(file['ti_{}/v1'.format(test_ind)])
    # for tind, t in enumerate(ts_LOP[1:], 1):
    for tind, t in enumerate(ts_LOP):
        # K = np.array(file['ti_{}/K'.format(tind)])
        # propagator = np.array(file['ti_{}/propagator'.format(tind)])  # F @ K @ Finv
        # E_test[tind] = np.sum(np.abs(propagator @ v1_test) ** 2.0)
        # KE_test[tind] = np.sum(np.abs(F_KE @ Finv @ propagator @ v1_test) ** 2.0)
        # ME_test[tind] = np.sum(np.abs(F_ME @ Finv @ propagator @ v1_test) ** 2.0)
        # U = np.array(file['ti_{}/U'.format(tind)])
        # Vh = np.array(file['ti_{}/Vh'.format(tind)])
        # s = np.array(file['ti_{}/s'.format(tind)])
        # Gs[tind] = s[0]**2.0
        LOP = np.array(file['ti_{}/LOP'.format(tind)])
        LOPs.append(LOP)

print('running t = {}'.format(ts_LOP[test_ind]))

z_basis = de.Chebyshev('z', Nz, interval=(-0.5 * Lz, 0.5 * Lz))
domain = de.Domain([z_basis], grid_dtype=np.complex128)
probvars = ['u','v','w','p','Ax','Ay','Az','Phi','omega_x','omega_y','Bx','By']

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
solver.stop_sim_time = ts_LOP[-1]
for i, var in enumerate(probvars):
    field = solver.state[var]
    field['c'] = LOPs[test_ind][int(i*len(LOPs[test_ind])/len(probvars)):int((i+1)*len(LOPs[test_ind])/len(probvars))]  # to generalize, 6 should become len(probvars)


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

profiles = solver.evaluator.add_file_handler(str(out_dir / 'profiles'), iter=10)
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

profiles_c = solver.evaluator.add_file_handler(str(out_dir / 'profiles_c'), iter=10)
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

dt = 2e-2
ts = np.arange(0.0, ts_LOP[-1] + dt, dt)
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

plt.semilogy(ts_LOP, Gs, c='k', label=r'$G(t_0)$')
plt.semilogy(ts_LOP, Gs[0]*np.exp(ts_LOP*2.0*np.max(np.real(eigenvalues[indx])))*Gs[-1]/np.exp(ts_LOP[-1]*2.0*np.max(np.real(eigenvalues[indx]))), '--', c='k', label=r'eigenvalue')
plt.semilogy(ts_LOP[test_ind], Gs[test_ind], 'x', c='k')
plt.semilogy(ts, energies, label=r'$E(t, t_0={:4.2f})$'.format(ts_LOP[test_ind]))
# plt.semilogy(ts_LOP, E_test, '--', label=r'$E$ via propagator')
plt.semilogy(ts, KEs, label=r'$KE(t, t_0={:4.2f})$'.format(ts_LOP[test_ind]))
# plt.semilogy(ts_LOP, KE_test, '--', label=r'$KE$ via propagator')
plt.semilogy(ts, MEs, label=r'$ME(t, t_0={:4.2f})$'.format(ts_LOP[test_ind]))
# plt.semilogy(ts_LOP, ME_test, '--', label=r'$ME$ via propagator')
plt.xlabel(r'$t, t_0$')
plt.ylabel(r'$G, E$')
plt.xlim((ts_LOP[0], ts_LOP[-1]))
plt.legend()
plt.savefig(plot_name+'_log.png')
plt.savefig(plot_name+'_log.pdf')
plt.yscale('linear')
plt.ylim(ymin=0)
plt.tight_layout()
plt.savefig(plot_name+'.png')
plt.savefig(plot_name+'.pdf')
plt.show()
