"""Solves the eigenvalue problem for KH instability of sinusoidal flow

Sets up the eigenvalue problem corresponding to the KH instability of a
sinusoidal flow profile in MHD with a uniform, flow-aligned magnetic 
field.

Can be used with viscosity/resistivity (ideal=False) or without 
(ideal=True). Currently only solves for 2D perturbations.

Note that all methods return frequencies or growth rates that are
normalized to the finger width/speed, like in Fig 3 of
Harrington & Garaud.

Methods
-------
Deln(k,n,delta)
    Calculates Delta_n from my LaTeX notes (for setting up KH EVP)
Lmat(delta, M2, Re, Rm, k, N, ideal=False)
    Constructs the linear operator whose eigenvalues are the
    complex frequencies for the KH eigenmodes, i.e., the modes
    are taken to go like f(t) ~ exp[i omega t], and the eigenvalues
    of this matrix are the various omega for different modes.
gamma_from_L(L)
    Calls numpy.linalg.eig to get the eigenvalues. Returns the
    growth rate gamma = -np.imag(omega) of the most unstable mode
omega_from_L(L)
    Returns the complex frequency instead (should just merge these
    two functions)
gamma_over_k(delta, M2, Re, Rm, ks, N, ideal=False)
    Calls gamfromL for each k in array ks, returns the resulting
    array of growth rates gamma[ks]
omega_over_k(...)
    Same as above but complex frequencies
"""
import numpy as np


def Deln(k, n, delta, finger_norm=False, k0=1.0):  # \Delta_n in my notes. So simple, probably shouldn't be a function
    if finger_norm:
        return k ** 2.0 + k0**2.0 * (n + delta)**2.0
    else:
        return k ** 2.0 + (n + delta) ** 2.0


def Lmat(delta, M2, Re, Rm, k, N, ideal=False):
    """Returns the linear operator for the KH instability

    Note that the eigenvalues are the complex frequencies, and the
    eigenvectors are the streamfunction and flux function mixed together,
    with the 2nth entry being the streamfunction at some wavenumber, and
    the 2n+1th being the flux function at a wavenumber.

    The eigenvalue problem this corresponds to is normalized to the flow speed and length scale, and background field.

    Parameters
    ----------
    delta : float
        This should be in the range 0 <= delta <= 0.5 and indicates
        the periodicity of the KH mode relative to the wavelength
        of the sinusoidal shear flow. See LaTeX notes; should
        probably be left at delta=0.0
    M2 : float
        This is 1/M_A^2 where M_A is the Alfven Mach number. 
        Equivalent to H_B^* in Harrington & Garaud.
    Re : float
        Reynolds number
    Rm : float
        Magnetic Reynolds number
    k : float
        Wavenumber in direction of flow
    N : int (ODD NUMBER)
        Numerical resolution in direction of shear. I've always 
        set it to an odd number. There's a few steps where I forgot to
        generalize to even or odd numbers for N.
    ideal : Bool, default=False
        Whether or not to set viscosity, resistivity -> 0
        (if True then Re and Rm don't matter)

    Returns
    -------
    L : 2N x 2N numpy array
        Matrix whose eigenvalues are complex frequencies of KH modes
    """
    diss = 1.0 - ideal  # =0 for ideal=True, =1 for ideal=False
    ns = list(range(-int((N - 1) / 2), int((N + 1) / 2), 1))
    ms = list(range(-N + 1, N + 1, 1))

    # the following few lines just sets up arrays of Delta_n
    delns = [Deln(k, n, delta) for n in ns]
    delns_m = np.zeros_like(ms, dtype=np.float64)
    for i, m in enumerate(ms):
        if m % 2 == 0:
            delns_m[i] = Deln(k, m / 2, delta)
        else:
            delns_m[i] = Deln(k, (m - 1) / 2, delta)

    M = 2 * N
    L = np.zeros((M, M), dtype=np.complex128)

    # first fill in the entries that aren't near the edges
    for i, m in enumerate(ms):
        deltan = delns_m[i]
        # deltanp1 = delns_m[i+2]
        # deltanm1 = delns_m[i-2]
        if i > 1 and i < len(ms) - 2:  # avoid entries near edges
            deltanp1 = delns_m[i + 2]
            deltanm1 = delns_m[i - 2]
            if m % 2 == 0:  # phi entries
                n = m / 2
                # phi_n, phi_n part
                L[i, i] = (1.0j) * (diss / Re) * deltan
                # phi_n, psi_n part
                L[i, i + 1] = M2 * k
                # phi_n, phi_n+1
                L[i, i + 2] = -k * (1 - deltanp1) / (2.0j * deltan)
                if not np.isfinite(L[i, i + 2]):
                    # Pretty sure I can get rid of this now -- I was debugging 0/0 errors, which I think only happen
                    # if you try to solve the system at k=0, which isn't interesting. And if it is, then the way to
                    # go about it is to multiply both sides of the linear system by Delta_n, and solve as a
                    # generalized eigenvalue problem
                    print(-k * (1 - deltanp1))
                    print(2.0j * deltan)
                # phi_n, phi_n-1
                L[i, i - 2] = k * (1 - deltanm1) / (2.0j * deltan)
            else:  # psi entries
                # psi_n, psi_n
                L[i, i] = (1.0j) * deltan * diss / Rm
                # psi_n, phi_n
                L[i, i - 1] = k
                # psi_n, psi_n+1
                L[i, i + 2] = k / (2.0j)
                # psi_n, psi_n-1
                L[i, i - 2] = -k / (2.0j)
    # now do the edges
    # first, the most negative phi
    L[0, 0] = (1.0j) * delns_m[0] * diss / Re
    L[0, 1] = M2 * k
    L[0, 2] = -k * (1 - delns_m[2]) / (2.0j * delns_m[0])
    # most negative psi
    L[1, 1] = (1.0j) * delns_m[1] * diss / Rm
    L[1, 0] = k
    L[1, 3] = k / (2.0j)
    # most positive phi
    L[-2, -2] = (1.0j) * delns_m[-2] * diss / Re
    L[-2, -1] = M2 * k
    L[-2, -4] = k * (1 - delns_m[-4]) / (2.0j * delns_m[-2])
    # most positive psi
    L[-1, -1] = (1.0j) * delns_m[-1] * diss / Rm
    L[-1, -2] = k
    L[-1, -3] = -k / (2.0j)
    return L


def gamma_from_L(L, withmode=False):
    w, v = np.linalg.eig(L)
    if withmode:
        ind = np.argmax(-np.imag(w))
        return [-np.imag(w[ind]), v[:, ind]]
    else:
        return np.max(-np.imag(w))


def omega_from_L(L):
    w, v = np.linalg.eig(L)
    wsort = w[np.argsort(-np.imag(w))]
    return wsort[-1]


def gamma_from_params(delta, M2, Re, Rm, k, N, ideal, withmode=False):
    L = Lmat(delta, M2, Re, Rm, k, N, ideal)
    return gamma_from_L(L, withmode)


def gamma_over_k(delta, M2, Re, Rm, ks, N, ideal=False):
    return [gamma_from_L(Lmat(delta, M2, Re, Rm, k, N, ideal)) for k in ks]


def omega_over_k(delta, M2, Re, Rm, ks, N, ideal=False):
    return [omega_from_L(Lmat(delta, M2, Re, Rm, k, N, ideal)) for k in ks]

