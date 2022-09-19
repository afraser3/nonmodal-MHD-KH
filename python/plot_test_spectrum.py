import h5py
from matplotlib import pyplot as plt
import sys
from pathlib import Path
plt.style.use('prl')

filename_2D = Path(sys.argv[-2])
filename_3D = Path(sys.argv[-1])

params = filename_2D.stem.split("_spectrum")[0]


df_2D = h5py.File(filename_2D, "r")
df_3D = h5py.File(filename_3D, "r")

plt.scatter(df_2D['spectrum'][:].real, df_2D['spectrum'][:].imag, label='2D', alpha=0.4)
plt.scatter(df_3D['spectrum'][:].real, df_3D['spectrum'][:].imag, label='2.5D', alpha=0.4)
plt.legend()
plt.xlim(-0.04,0)
plt.xlabel('real')
plt.ylabel('imag')
plt.tight_layout()

plt.savefig(f'spectrum_2D_2.5D_{params}.png', dpi=300)
