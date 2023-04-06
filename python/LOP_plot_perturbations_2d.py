# Reads the files made by LOP_save_mats.py and plots the perturbations
import numpy as np
from matplotlib import pyplot as plt
import h5py
from dedalus import public as de

Nz = 512
Reynolds = 250.0
Pm = 1.0
mReynolds = Pm * Reynolds
Lz = 10.0*np.pi
kx = 0.2
MA = 1.25
MA2 = MA**2.0
k = 100
dense = True
reject_spurious = True
if dense:
    fname = "LOPs/Re{}/MA{}/Pm{}/Nz{}_dense.h5".format(Reynolds, MA, Pm, Nz)
else:
    fname = "LOPs/Re{}/MA{}/Pm{}/Nz{}_k{}.h5".format(Reynolds, MA, Pm, Nz, k)

with h5py.File(fname, mode='r') as file:
    K = np.array(file['ti_10/K'])
    U = np.array(file['ti_10/U'])
    Vh = np.array(file['ti_10/Vh'])
    s = np.array(file['ti_10/s'])
    eigenvalues = np.array(file['setup_mats/eigenvalues'])
    indx = np.argmax(np.real(eigenvalues))
    evec = np.array(file['setup_mats/eigenvectors'][:, indx])
    F = np.array(file['setup_mats/F'])
    R = np.array(file['setup_mats/R'])
    Q = np.array(file['setup_mats/Q'])
    pre_right = np.array(file['setup_mats/pre_right'])
    print(np.array(file['scales/ts']))
    Nz = np.array(file['scales/Nz'])
    Lz = np.array(file['scales/Lz'])
    kx = np.array(file['scales/kx'])

v1 = Vh.conj()[0]  # the optimal perturbation in the F'ed and Q'ed basis
v1F = np.linalg.inv(F) @ v1  # change of basis away from the F'ed basis
# now I want to say that R is the change-of-basis matrix between the Q and eigenvector bases, and so I should be doing
# R^-1, but that doesn't give the right dimensions.
# I'm pretty sure that v1F is v1 in the basis of the columns of Q, so I should be able to do Q @ v1F to
# get the perturbation in terms of Chebyshev expansion coefficients
LOP = Q @ v1F
LOP = pre_right @ LOP

evec_phase = evec[0]/np.abs(evec[0])
LOP_phase = LOP[0]/np.abs(LOP[0])
evec = evec/evec_phase
LOP = LOP/LOP_phase

z_basis = de.Chebyshev('z', Nz, interval=(-0.5 * Lz, 0.5 * Lz))
x_basis = de.Fourier('x', 32, interval=(0.0, 2.0 * np.pi / kx))
domain = de.Domain([x_basis, z_basis], grid_dtype=np.float64)
# domain = de.Domain([z_basis], grid_dtype=np.complex128)
x = domain.grid(0)[:, 0]
z = domain.grid(1)[0]
# z = domain.grid(0)
phi_LOP = domain.new_field(name='phi_LOP')
psi_LOP = domain.new_field(name='psi_LOP')
phi_evec = domain.new_field(name='phi_evec')
psi_evec = domain.new_field(name='psi_evec')
# phi_z = domain.new_field(name='phi_z')
#phi_LOP['c'][1, :] = LOP[::6]
#psi_LOP['c'][1, :] = LOP[2::6]
phi_LOP['c'][1, :] = LOP[:int(len(LOP)/6)]
psi_LOP['c'][1, :] = LOP[int(2*len(LOP)/6):int(3*len(LOP)/6)]
phi_evec['c'][1, :] = evec[:int(len(evec)/6)]  # evec[::6]
psi_evec['c'][1, :] = evec[int(2*len(evec)/6):int(3*len(evec)/6)]  # [2::6]
# phi['c'] = LOP[::6]
phi_LOP_out = phi_LOP['g']
psi_LOP_out = psi_LOP['g']
phi_evec_out = phi_evec['g']
psi_evec_out = psi_evec['g']
# phi_z['c'] = LOP[1::6]
# phi_z_out = phi_z['g']

plt.subplot(2, 2, 1)
plt.contourf(x, z, phi_LOP_out.T)
plt.title(r'LOP $\phi$')
#plt.ylim((-7.5, 7.5))
plt.colorbar()
# plt.plot(z, np.real(phi_out))
# plt.plot(z, np.real(phi_z_out))

plt.subplot(2, 2, 2)
plt.contourf(x, z, psi_LOP_out.T)
plt.title(r'LOP $\psi$')
#plt.ylim((-7.5, 7.5))
plt.colorbar()

plt.subplot(2, 2, 3)
plt.contourf(x, z, phi_evec_out.T)
plt.title(r'Eigenmode $\phi$')
#plt.ylim((-7.5, 7.5))
plt.colorbar()

plt.subplot(2, 2, 4)
plt.contourf(x, z, psi_evec_out.T)
plt.title(r'Eigenmode $\psi$')
#plt.ylim((-7.5, 7.5))
plt.colorbar()
plt.tight_layout()
plt.show()
#plt.savefig('LOP_eigenmode_comparison_MA1.25.png')
