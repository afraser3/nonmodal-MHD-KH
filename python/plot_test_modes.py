"""

mode is mode number for 2D
mode_3D is mode for 3D.

By default, test_implementation_2D.py returns 14 modes; test_implementation_3D.py returns 4.

The first two in 3D are the least stable ones present in 2D; the second two are the 2.5D modes with no counterpart in 2D.

"""

import h5py
from matplotlib import pyplot as plt
import sys
from pathlib import Path
plt.style.use('prl')

filename_2D = Path(sys.argv[-2])
filename_3D = Path(sys.argv[-1])

params = filename_2D.stem.split("_spectrum")[0]

mode = 13
mode_3D = 0

df_2D = h5py.File(filename_2D, "r")
df_3D = h5py.File(filename_3D, "r")

z_2D = df_2D['scales/z'][:]
z_3D = df_3D['scales/z'][:]

eval_2D = df_2D['spectrum'][mode]
eval_3D = df_3D['spectrum'][mode_3D]
print(df_3D['spectrum'][:])
print(df_2D['spectrum'][:])

u_2D = df_2D[f'/mode_{mode:d}/u'][:]
w_2D = df_2D[f'/mode_{mode:d}/w'][:]
Bx_2D = df_2D[f'/mode_{mode:d}/Bx'][:]
Bz_2D = df_2D[f'/mode_{mode:d}/Bz'][:]

ampl_2D = u_2D.conj()
ampl_2D /= (ampl_2D*u_2D).real.max()

u_3D = df_3D[f'/mode_{mode_3D:d}/u'][:]
v_3D = df_3D[f'/mode_{mode_3D:d}/v'][:]
w_3D = df_3D[f'/mode_{mode_3D:d}/w'][:]
Bx_3D = df_3D[f'/mode_{mode_3D:d}/Bx'][:]
By_3D = df_3D[f'/mode_{mode_3D:d}/By'][:]
Bz_3D = df_3D[f'/mode_{mode_3D:d}/Bz'][:]

ampl_3D = u_3D.conj()
ampl_3D /= (ampl_3D*u_3D).real.max()
fig = plt.figure(figsize=(8,5))

gs = fig.add_gridspec(2, 3, hspace=0.1, wspace=0)
(u_ax, v_ax, w_ax), (Bx_ax, By_ax, Bz_ax) = gs.subplots(sharex='col', sharey='row')
u_ax.plot(z_2D, (ampl_2D*u_2D).real, label=f'2D real',)
u_ax.plot(z_3D, (ampl_3D*u_3D).real, '-.', label=f'3D real',)
u_ax.set_prop_cycle(None)
u_ax.plot(z_2D, (ampl_2D*u_2D).imag, '--', label='imag')
u_ax.plot(z_3D, (ampl_3D*u_3D).imag, ':', label='imag')
u_ax.set_ylabel("u")


v_ax.plot([],[], label='2D real')
v_ax.plot(z_3D, (ampl_3D*v_3D).real, '-.', label='3D real')
v_ax.set_prop_cycle(None)
v_ax.plot([],[], '--', label='imag')
v_ax.plot([],[], ':', label='imag')
v_ax.legend(fontsize=10,ncol=2)
w_ax.plot(z_2D, (ampl_2D*w_2D).real,)
w_ax.plot(z_3D, (ampl_3D*w_3D).real,'-.',)
w_ax.set_prop_cycle(None)
w_ax.plot(z_2D, (ampl_2D*w_2D).imag,'--')
w_ax.plot(z_3D, (ampl_3D*w_3D).imag,':')

u_ax.set_title("x")
v_ax.set_title("y")
w_ax.set_title("z")

Bx_ax.plot(z_2D, (ampl_2D*Bx_2D).real)
Bx_ax.plot(z_3D, (ampl_3D*Bx_3D).real,'-.')
Bx_ax.set_prop_cycle(None)
Bx_ax.plot(z_2D, (ampl_2D*Bx_2D).imag,'--')
Bx_ax.plot(z_3D, (ampl_3D*Bx_3D).imag,'--')

Bx_ax.set_ylabel("B")
for ax in [Bx_ax, By_ax, Bz_ax]:
    ax.set_xlabel("z")
By_ax.plot([],[])
By_ax.plot(z_3D, (ampl_3D*By_3D).real, '-.')


Bz_ax.plot(z_2D, (ampl_2D*Bz_2D).real)
Bz_ax.plot(z_3D, (ampl_3D*Bz_3D).real, '-.')
Bz_ax.set_prop_cycle(None)
Bz_ax.plot(z_2D, (ampl_2D*Bx_2D).imag,'--')
Bz_ax.plot(z_3D, (ampl_3D*Bx_3D).imag,':')


fig.savefig(f'emodes_mode_{mode:d}_{params}.pdf')
