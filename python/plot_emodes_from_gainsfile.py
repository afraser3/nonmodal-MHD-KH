"""
The plot_modes stuff in scan_linear_optimals_3d.py is helpful but it creates a bunch of different pdf's that can be
tedious to sift through. The goal here is to create, in each kx/ky folder, a single multi-page pdf where each page
is a plot of a mode, including as a subplot the spectra so you can see where in the spectrum this mode lies
"""
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import h5py
from dedalus import public as de

basedir = 'LOPs_3D/run_kx-log-scan-saveLOPs-saveEvecs/kxi_0015/kyi_0000/'
Legendre = False  # should fix this to be read from the config file or something
gainsfile = basedir+'gains_mats.h5'
with h5py.File(gainsfile, mode='r') as file:
    kx = np.array(file['scales/kx'])
    ky = np.array(file['scales/ky'])
    Nz = np.array(file['scales/Nz'])
    Lz = np.array(file['scales/Lz'])
    z = np.array(file['scales/z'])
    evals = np.array(file['setup_mats/eigenvalues'])
    evecs = np.array(file['setup_mats/eigenvectors'])

with PdfPages(basedir + 'evecs.pdf') as pdf:
    if Legendre:
        z_basis = de.Legendre('z', Nz, interval=(-0.5 * Lz, 0.5 * Lz))
    else:
        z_basis = de.Chebyshev('z', Nz, interval=(-0.5 * Lz, 0.5 * Lz))
    x_basis = de.Fourier('x', 32, interval=(0.0, 2.0 * np.pi / kx))
    y_basis = de.Fourier('x', 32, interval=(0.0, 2.0 * np.pi / ky))
    domain = de.Domain([x_basis, y_basis, z_basis], grid_dtype=np.float64)
    x = domain.grid(0)[:, 0, 0]
    y = domain.grid(1)[0, :, 0]
    z = domain.grid(2)[0, 0, :]
    probvars = ['u', 'v', 'w', 'p', 'Ax', 'Ay', 'Az', 'Phi', 'omega_x', 'omega_y', 'Bx', 'By']
    evec_fields = [domain.new_field(name=probvar) for probvar in probvars]
    Bz = domain.new_field(name='Bz')
    for evi, eval in enumerate(evals):
        if evi%10==0:
            print('evi = ', evi)
        evec = evecs[:, evi]
        evec_phase = evec[0] / np.abs(evec[0])
        evec = evec / evec_phase
        for ind, field in enumerate(evec_fields):
            # select the part of the array that corresponds to each field
            slice_start = int(ind*len(evec)/len(probvars))
            slice_stop = int((ind+1)*len(evec)/len(probvars))
            field['c'][1, 1, :] = evec[slice_start:slice_stop]
        Bz['c'] = 1j*kx*evec_fields[5]['c'] - 1j*ky*evec_fields[4]['c']  # Bz = dx(Ay) - dy(Ax)

        yind = 0  # where in y to slice)
        scale = 0.8
        fig = plt.figure(figsize=(scale*16, scale*10))
        plt.subplot(2, 4, 1)
        plt.pcolormesh(x, z, evec_fields[0]['g'][:, yind, :].T, rasterized=True)
        # plt.title(r'$u$, $t={}$'.format(LOP_ts[tind]))
        plt.title(r'$u$')
        plt.ylabel(r'$z$')
        plt.colorbar()

        plt.subplot(2, 4, 2)
        plt.pcolormesh(x, z, evec_fields[1]['g'][:, yind, :].T, rasterized=True)
        plt.title(r'$v$')
        # plt.ylabel(r'$z$')
        plt.colorbar()

        plt.subplot(2, 4, 3)
        plt.pcolormesh(x, z, evec_fields[2]['g'][:, yind, :].T, rasterized=True)
        plt.title(r'$w$')
        # plt.ylabel(r'$z$')
        plt.colorbar()

        plt.subplot(2, 4, 5)
        plt.pcolormesh(x, z, evec_fields[10]['g'][:, yind, :].T, rasterized=True)
        plt.title(r'$B_x$')
        plt.ylabel(r'$z$')
        plt.xlabel(r'$x$')
        plt.colorbar()

        plt.subplot(2, 4, 6)
        plt.pcolormesh(x, z, evec_fields[11]['g'][:, yind, :].T, rasterized=True)
        plt.title(r'$B_y$')
        # plt.ylabel(r'$z$')
        plt.xlabel(r'$x$')
        plt.colorbar()

        plt.subplot(2, 4, 7)
        plt.pcolormesh(x, z, Bz['g'][:, yind, :].T, rasterized=True)
        plt.title(r'$B_z$')
        # plt.ylabel(r'$z$')
        plt.xlabel(r'$x$')
        plt.colorbar()

        plt.subplot(2, 4, 4)
        plt.plot(np.real(evals), np.imag(evals), '.')
        plt.plot(np.real(eval), np.imag(eval), 'x')
        plt.xlabel(r'growth rate, $Re(\gamma)$')
        plt.ylabel(r'frequency, $Im(\gamma)$')

        plt.suptitle(r'$k_x={:1.2e}$, $\gamma={:1.3e}$'.format(kx, eval))

        plt.tight_layout()
        # plt.show()
        # plt.savefig(str(kydir / 'LOP_tind{}.png'.format(tind)))
        # plt.savefig(str(kydir / 'LOP_tind{}.pdf'.format(tind)))
        pdf.savefig()
        plt.close()

