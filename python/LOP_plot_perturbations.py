# Reads the files made by LOP_save_mats.py and plots the perturbations
import numpy as np
from matplotlib import pyplot as plt
import h5py
from dedalus import public as de

Nz = 256
Reynolds = 250.0
Pm = 1.0
mReynolds = Pm * Reynolds
Lz = 10.0*np.pi
kx = 0.2
MA = 0.75  # For this value, the system is stable.
MA2 = MA**2.0
k = 100
dense = False
reject_spurious = True
if dense:
    fname = "LOPs/Re{}/MA{}/Pm{}/Nz{}_dense.h5".format(Reynolds, MA, Pm, Nz)
else:
    fname = "LOPs/Re{}/MA{}/Pm{}/Nz{}_k{}.h5".format(Reynolds, MA, Pm, Nz, k)

with h5py.File(fname, mode='r') as file:
    K = np.array(file['ti_2/K'])
    U = np.array(file['ti_2/U'])
    Vh = np.array(file['ti_2/Vh'])
    s = np.array(file['ti_2/s'])
    F = np.array(file['setup_mats/F'])
    R = np.array(file['setup_mats/R'])
    Q = np.array(file['setup_mats/Q'])
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

z_basis = de.Chebyshev('z', Nz, interval=(-0.5 * Lz, 0.5 * Lz))
x_basis = de.Fourier('x', 32, interval=(0.0, 2.0 * np.pi / kx))
domain = de.Domain([x_basis, z_basis], grid_dtype=np.float64)
# domain = de.Domain([z_basis], grid_dtype=np.complex128)
x = domain.grid(0)[:, 0]
z = domain.grid(1)[0]
# z = domain.grid(0)
phi = domain.new_field(name='phi')
# phi_z = domain.new_field(name='phi_z')
phi['c'][1, :] = LOP[::6]
# phi['c'] = LOP[::6]
phi_out = phi['g']
# phi_z['c'] = LOP[1::6]
# phi_z_out = phi_z['g']

plt.contourf(x, z, phi_out.T)
plt.colorbar()
# plt.plot(z, np.real(phi_out))
# plt.plot(z, np.real(phi_z_out))
plt.show()
