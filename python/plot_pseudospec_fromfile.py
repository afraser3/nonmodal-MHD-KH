"""
Copying this from an old script I had lying around from grad school.
I did very little testing of the viscous/resistive system, so it needs to be checked.
(I have verified that it returns a reasonable eigenvalue for the current set of parameters.)
It could definitely also use cleaning up.
"""
import h5py
import numpy as np
from matplotlib import pyplot as plt

Reynolds = 1000.0
mReynolds = 500.0
kx = 0.2
MA = 1.2

Pm = mReynolds/Reynolds
Nz = 512
Lz = 10.0*np.pi

k = 800

filepath = 'saved_spectra/Pm{}/Re{}/MA{}/kx{}/'.format(Pm, Reynolds, MA, kx)
filename = 'pseudospec_Nz{}_Lz{}pi_k{}.h5'.format(Nz, Lz/np.pi, k)
plotname = 'pseudospec_Nz{}_Lz{}pi_k{}.pdf'.format(Nz, Lz/np.pi, k)

with h5py.File(filepath+filename, 'r') as file:
    evalues = np.array(file['evalues/evalues'])
    ps_real = np.array(file['pseudospec/ps_real'])
    ps_imag = np.array(file['pseudospec/ps_imag'])
    pseudospectrum = np.array(file['pseudospec/pseudospectrum'])
plt.plot(np.real(evalues), np.imag(evalues), '.', c='k')
# plt.legend()
plt.contour(ps_real, ps_imag, np.log10(pseudospectrum), levels=np.arange(-8, 0))
plt.colorbar(label=r'$\log_{10} (\epsilon)$')
plt.ylim((ps_imag[0], ps_imag[-1]))
plt.xlim((ps_real[0], ps_real[-1]))
plt.axhline(0, color='k', alpha=0.2)
plt.xlabel('real (oscillating) frequency')
plt.ylabel('growth rate')
plt.title(r'(kx, MA, Re, Pm, k) = ({}, {}, {}, {}, {})'.format(kx, MA, Reynolds, Pm, k))
plt.savefig(filepath+plotname)
# plt.show()
