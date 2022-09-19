# Reads the files made by LOP_save_mats.py and plots the perturbations
import numpy as np
from matplotlib import pyplot as plt
import h5py
from dedalus import public as de
import logging
logger = logging.getLogger(__name__)

test_ind = 20

Nz = 256
Reynolds = 250.0
Pm = 1.0
mReynolds = Pm * Reynolds
Lz = 10.0*np.pi
kx = 0.2
MA = 0.85
MA2 = MA**2.0
k = 100
dense = True
reject_spurious = True
if dense:
    fname = "LOPs/Re{}/MA{}/Pm{}/Nz{}_dense.h5".format(Reynolds, MA, Pm, Nz)
    out_dir = "LOPs/Re{}/MA{}/Pm{}/Nz{}_dense/".format(Reynolds, MA, Pm, Nz)
else:
    fname = "LOPs/Re{}/MA{}/Pm{}/Nz{}_k{}.h5".format(Reynolds, MA, Pm, Nz, k)
    out_dir = "LOPs/Re{}/MA{}/Pm{}/Nz{}_k{}/".format(Reynolds, MA, Pm, Nz, k)

with h5py.File(fname, mode='r') as file:
    eigenvalues = np.array(file['setup_mats/eigenvalues'])
    indx = np.argmax(np.real(eigenvalues))
    evec = np.array(file['setup_mats/eigenvectors'][:, indx])
    evec_phase = evec[0] / np.abs(evec[0])
    evec = evec / evec_phase
    F = np.array(file['setup_mats/F'])
    Finv = np.linalg.inv(F)
    F_KE = np.array(file['setup_mats/F_KE'])
    F_ME = np.array(file['setup_mats/F_ME'])
    R = np.array(file['setup_mats/R'])
    Q = np.array(file['setup_mats/Q'])
    pre_right = np.array(file['setup_mats/pre_right'])
    ts_LOP = np.array(file['scales/ts'])
    Nz = np.array(file['scales/Nz'])[()]
    Lz = np.array(file['scales/Lz'])[()]
    kx = np.array(file['scales/kx'])[()]
    LOPs = []
    Gs = np.zeros_like(ts_LOP)
    E_test = np.zeros_like(ts_LOP)
    KE_test = np.zeros_like(ts_LOP)
    ME_test = np.zeros_like(ts_LOP)
    v1_test = np.array(file['ti_{}/Vh'.format(test_ind)]).conj()[0]
    # for tind, t in enumerate(ts_LOP[1:], 1):
    for tind, t in enumerate(ts_LOP):
        # K = np.array(file['ti_{}/K'.format(tind)])
        propagator = np.array(file['ti_{}/propagator'.format(tind)])  # F @ K @ Finv
        E_test[tind] = np.sum(np.abs(propagator @ v1_test) ** 2.0)
        KE_test[tind] = np.sum(np.abs(F_KE @ Finv @ propagator @ v1_test) ** 2.0)
        ME_test[tind] = np.sum(np.abs(F_ME @ Finv @ propagator @ v1_test) ** 2.0)
        # U = np.array(file['ti_{}/U'.format(tind)])
        Vh = np.array(file['ti_{}/Vh'.format(tind)])
        s = np.array(file['ti_{}/s'.format(tind)])
        Gs[tind] = s[0]**2.0
        v1 = Vh.conj()[0]  # the optimal perturbation in the F'ed and Q'ed basis
        v1F = Finv @ v1  # change of basis away from the F'ed basis
        # now I want to say that R is the change-of-basis matrix between the Q and eigenvector bases, and so I should
        # be doing R^-1, but that doesn't give the right dimensions.
        # I'm pretty sure that v1F is v1 in the basis of the columns of Q, so I should be able to do Q @ v1F to
        # get the perturbation in terms of Chebyshev expansion coefficients
        LOP = Q @ v1F
        LOP = pre_right @ LOP
        LOP_phase = LOP[0]/np.abs(LOP[0])
        LOP = LOP/LOP_phase
        LOPs.append(LOP)

z_basis = de.Chebyshev('z', Nz, interval=(-0.5 * Lz, 0.5 * Lz))
domain = de.Domain([z_basis], grid_dtype=np.complex128)
z = domain.grid(0)
probvars = ['phi', 'phi_z', 'psi', 'psi_z', 'uz', 'uzz']
problem = de.IVP(domain, variables=probvars, time='t')

a = 1.0
nccU = domain.new_field(name='U')
nccU['g'] = np.tanh(z / a)
problem.parameters['U'] = nccU
nccDU = domain.new_field(name='DU')
nccDU['g'] = 1 / (a * np.cosh(z / a) ** 2)
problem.parameters['DU'] = nccDU

problem.parameters['MA2'] = MA2  # Alfven Mach number squared
problem.parameters['kx'] = kx
problem.parameters['Lz'] = Lz
problem.parameters['Re'] = Reynolds
problem.parameters['Rm'] = mReynolds
problem.substitutions['dx(A)'] = "1.0j*kx*A"
problem.substitutions['u'] = 'phi_z'
problem.substitutions['w'] = "-dx(phi)"
problem.substitutions['Bz'] = "-dx(psi)"
problem.substitutions['Bx'] = "psi_z"
problem.substitutions['Jz'] = "dx(dx(psi)) + dz(Bx)"
problem.substitutions['zeta'] = "dx(dx(phi)) + uz"
problem.substitutions['zetaz'] = "dx(dx(u)) + uzz"
problem.add_equation("dt(zeta) - 1/Re*(dx(dx(zeta)) + dz(zetaz)) + U*dx(zeta) - dx(phi)*dz(DU) - 1/MA2*(dx(Jz)) = 0")
# problem.add_equation("dt(zeta) - 1/Re*(dx(dx(zeta)) + dz(zetaz)) - 1/MA2*(dx(Jz)) = -U*dx(zeta) + dx(phi)*dz(DU)")
problem.add_equation("uz - dz(u) = 0")
problem.add_equation("uzz - dz(uz) = 0")
problem.add_equation("dt(psi) - 1/Rm*Jz + U*dx(psi) - dx(phi) = 0")
# problem.add_equation("dt(psi) - 1/Rm*Jz - dx(phi) = -U*dx(psi)")
problem.add_equation("phi_z - dz(phi) = 0")
problem.add_equation("psi_z - dz(psi) = 0")

problem.add_bc("left(w) = 0")
problem.add_bc("right(w) = 0")
problem.add_bc("left(dz(w)) = 0")  # no-slip, co-moving (with base flow U_0) boundaries
problem.add_bc("right(dz(w)) = 0")  # (does free-slip make more sense?)
problem.add_bc("right(Bz) = 0")  # perfectly conducting boundaries
problem.add_bc("left(Bz) = 0")

logger.info("done with BCs")

solver = problem.build_solver(de.timesteppers.RK443)
solver.stop_sim_time = ts_LOP[-1]
for i, var in enumerate(probvars):
    field = solver.state[var]
    field['c'] = LOPs[test_ind][int(i*len(LOPs[test_ind])/6):int((i+1)*len(LOPs[test_ind])/6)]  # to generalize, 6 should become len(probvars)


def energy_norm(X1):
    u1 = X1['phi_z']
    w1 = -1.0j * kx * X1['phi']
    bx1 = X1['psi_z']
    bz1 = -1.0j * kx * X1['psi']
    field = (np.abs(u1) ** 2.0 + np.abs(w1) ** 2.0 + (1.0 / MA2) * (
                np.abs(bx1) ** 2.0 + np.abs(bz1) ** 2.0)).evaluate().integrate()
    return field['g'][0]


def kinetic_energy(X1):
    u1 = X1['phi_z']
    w1 = -1.0j * kx * X1['phi']
    field = (np.abs(u1) ** 2.0 + np.abs(w1) ** 2.0).evaluate().integrate()
    return field['g'][0]


def magnetic_energy(X1):
    bx1 = X1['psi_z']
    bz1 = -1.0j * kx * X1['psi']
    field = ((1.0 / MA2) * (np.abs(bx1) ** 2.0 + np.abs(bz1) ** 2.0)).evaluate().integrate()
    return field['g'][0]


dt = 1e-2
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
plt.semilogy(ts_LOP[test_ind], Gs[test_ind], 'x', c='k')
plt.semilogy(ts, energies, label=r'$E(t, t_0={:4.2f})$ via IVP'.format(ts_LOP[test_ind]))
plt.semilogy(ts_LOP, E_test, '--', label=r'$E$ via propagator')
plt.semilogy(ts, KEs, label=r'$KE(t, t_0={:4.2f})$'.format(ts_LOP[test_ind]))
plt.semilogy(ts_LOP, KE_test, '--', label=r'$KE$ via propagator')
plt.semilogy(ts, MEs, label=r'$ME(t, t_0={:4.2f})$'.format(ts_LOP[test_ind]))
plt.semilogy(ts_LOP, ME_test, '--', label=r'$ME$ via propagator')
plt.xlabel(r'$t, t_0$')
plt.ylabel(r'$G, E$')
plt.xlim((ts_LOP[0], ts_LOP[-1]))
plt.legend()
plt.tight_layout()
plt.yscale('linear')
plt.savefig('LOP_IVP_test3.pdf')
# plt.show()
