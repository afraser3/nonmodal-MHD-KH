"""
Copying this from an old script I had lying around from grad school.
I did very little testing of the viscous/resistive system, so it needs to be checked.
(I have verified that it returns a reasonable eigenvalue for the current set of parameters.)
It could definitely also use cleaning up.
"""
# from eigentools import Eigenproblem
from dedalus import public as de
import numpy as np

Reynolds = 500.0
mReynolds = 500.0
Nz = 128
Lz = 10.0*np.pi
kx = 0.2
MA = 10.0

solvedense = False
nev = 1
findstable = False
energy_norm = True
savemodes = False
basefolder = 'savedmodes/'


def compute_energies(solver):  # Took this from somewhere on Keaton Burns' GitHub
    # Construct energy operator
    Lz = solver.problem.parameters['Lz']
    kx = solver.problem.parameters['kx']
    Lx = 2.0 * np.pi / kx
    MA2 = solver.problem.parameters['MA2']
    phi = solver.state['phi']
    psi = solver.state['psi']

    phi_z_op = de.operators.differentiate(phi, 'z')
    psi_z_op = de.operators.differentiate(psi, 'z')

    E_op = (Lx / 2.0) * de.operators.integrate(
        np.abs(kx ** 2 * phi * phi) + np.abs(phi_z_op * phi_z_op) + (1.0 / MA2) * (
                    np.abs(kx ** 2 * psi * psi) + np.abs(psi_z_op * psi_z_op)), 'z')
    # Evaluate energy for each mode
    N = len(solver.eigenvalues)
    energies = np.zeros(N)
    for i in range(N):
        solver.set_state(i)
        energies[i] = np.abs(E_op.evaluate()['c'][0])
    return energies


# def KHevp(Nz, Lz, kx, MA, Reynolds=np.inf, mReynolds=np.inf, bounds=0, splitdomain=False, plotmodes=False,
          # solvedense=False, nev=1, doleft=False, findstable=False, test=False, energy_norm=True, basefolder=''):

if Reynolds == np.inf:
    inviscid = True
else:
    inviscid = False
if mReynolds == np.inf:
    ideal = True
else:
    ideal = False
z_basis = de.Chebyshev('z', Nz, interval=(-0.5 * Lz, 0.5 * Lz))
domain = de.Domain([z_basis], grid_dtype=np.complex128)
probvars = ['phi', 'phi_z', 'psi', 'psi_z']
if inviscid == False:
    probvars.append('uz')
    probvars.append('uzz')
else:
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

problem.parameters['MA2'] = MA ** 2.0  # Alfven Mach number squared
problem.parameters['kx'] = kx
problem.parameters['Lz'] = Lz
if inviscid == False:
    problem.parameters['Re'] = Reynolds
if ideal == False:
    problem.parameters['Rm'] = mReynolds
problem.substitutions['dx(A)'] = "1.0j*kx*A"
problem.substitutions['dt(A)'] = "1.0j*omega*A"
problem.substitutions['u'] = 'phi_z'
problem.substitutions['w'] = "-dx(phi)"
problem.substitutions['Bz'] = "-dx(psi)"
problem.substitutions['Bx'] = "psi_z"
problem.substitutions['Jz'] = "dx(dx(psi)) + dz(Bx)"
if inviscid == False:
    problem.substitutions['zeta'] = "dx(dx(phi)) + uz"
    problem.substitutions['zetaz'] = "dx(dx(u)) + uzz"
if inviscid == False:
    problem.add_equation(
        "dt(zeta) - 1/Re*(dx(dx(zeta)) + dz(zetaz)) + U*dx(zeta) - dx(phi)*dz(DU) - 1/MA2*(dx(Jz)) = 0")
    problem.add_equation("uz - dz(u) = 0")
    problem.add_equation("uzz - dz(uz) = 0")
else:
    problem.add_equation("dt(zeta) + U*dx(zeta) - dx(phi)*dz(DU) - 1/MA2*(dx(Jz)) = 0")
    problem.add_equation("zeta - dx(dx(phi)) - dz(phi_z) = 0")
if ideal == False:
    problem.add_equation("dt(psi) - 1/Rm*Jz + U*dx(psi) - dx(phi) = 0")
else:
    problem.add_equation("dt(psi) + U*dx(psi) - dx(phi) = 0")
problem.add_equation("phi_z - dz(phi) = 0")
problem.add_equation("psi_z - dz(psi) = 0")

problem.add_bc("left(w) = 0")
problem.add_bc("right(w) = 0")
if inviscid == False:
    problem.add_bc("left(dz(w)) = 0")
    problem.add_bc("right(dz(w)) = 0")
problem.add_bc("right(Bz) = 0")
problem.add_bc("left(Bz) = 0")

print("done with BCs")

solver = problem.build_solver()
if solvedense == True:
    solver.solve_dense(solver.pencils[0])
else:  # very rough guess
    if kx < 0.5:
        guess = -1.0j * 0.5 * kx  # * np.sqrt(1.0 - (2.0/MA)**2)  # from palotti eq 7
    else:
        guess = -1.0j * 0.5 * (1.0 - kx)  # * np.sqrt(1.0-(2.0/MA)**2)

    if findstable == True:
        guess = -guess
    solver.solve_sparse(solver.pencils[0], nev, 2.0 * guess, rebuild_coeffs=False)

ev = solver.eigenvalues
vecs = solver.eigenvectors

# normalize phase of eigenvectors
phase = lambda Z: Z / np.abs(Z)
solver.eigenvectors /= phase(solver.eigenvectors[0:1, :])
if energy_norm == True:
    energies = compute_energies(solver)
    solver.eigenvectors /= np.sqrt(energies)
    vecs = solver.eigenvectors

# end_time = time.time()
# logger.info('Run time: %.2f sec' % (end_time - start_time))


if savemodes == True:
    sort = True
    if sort == True:
        # Remove nan and inf eigenvalues and sort them w/ largest gamma first
        # -inf first, then smallest gamma, nans/+infs last
        # order = np.argsort(-np.imag(ev))
        # sort instead by increasing |omega|
        order = np.argsort(np.abs(np.real(ev)))
        while np.isfinite(ev[order[0]]) == False:  # take leading -infs
            order = np.append(order[1:], order[0])  # and put them in the back
        # Finally, reorganize eigenvalues and vectors according to this order
        # ev = ev[order]
        # vecs = vecs[order]
        # Trim trailing nans and infs
        ##vecs = vecs[np.isfinite(ev)]
        ##ev = ev[np.isfinite(ev)]
        # order = order[np.isfinite(ev[order])]

    if solvedense == False and findstable == True:
        analysisc = solver.evaluator.add_file_handler(basefolder + 'Eigenmodes_c_stable', iter=np.inf)
        analysisc.add_system(solver.state, layout='c')

        analysisg = solver.evaluator.add_file_handler(basefolder + 'Eigenmodes_g_stable', iter=np.inf)
        analysisg.add_system(solver.state, layout='g')

        for i in range(len(order)):
            # if np.isfinite(ev[i]):
            solver.set_state(order[i])
            if i == 0:
                analysisc.create_current_file()
                analysisg.create_current_file()
            solver.evaluator.evaluate_handlers([analysisc, analysisg], world_time=kx,
                                               wall_time=np.imag(ev[order[i]]), sim_time=np.real(ev[order[i]]),
                                               timestep=order[i], iteration=i)

    else:
        analysisc = solver.evaluator.add_file_handler(basefolder + 'Eigenmodes_c', iter=np.inf)
        analysisc.add_system(solver.state, layout='c')

        analysisg = solver.evaluator.add_file_handler(basefolder + 'Eigenmodes_g', iter=np.inf)
        analysisg.add_system(solver.state, layout='g')

        for i in range(len(order)):
            # if np.isfinite(ev[i]):
            solver.set_state(order[i])
            if i == 0:
                analysisc.create_current_file()
                analysisg.create_current_file()
            solver.evaluator.evaluate_handlers([analysisc, analysisg], world_time=kx,
                                               wall_time=np.imag(ev[order[i]]), sim_time=np.real(ev[order[i]]),
                                               timestep=order[i], iteration=i)

