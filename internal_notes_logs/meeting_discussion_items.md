## June 3rd, 2021
### Progress
- Calculated my first pseudospectrum! The functioning script is in nonmodal-MHD-KH/python/MHD_OrrSomm_pseudospec.py
  
### Questions/issues/discussion
- Can we talk about how Eigentools calculates pseudospectra under the hood?
  - Noticed calc_ps calls EP.solve(sparse=True). If I first call EP.solve(sparse=False), check the % difference between 
    evalues_low and _high for the unstable mode, then check it again after calling calc_ps, the % difference goes up,
    causing the mode to go from not-rejected to rejected.
- Dissipative or dissipationless spectra?
- Do you usually calculate pseudospectra at all unstable wavenumbers? Most-unstable wavenumber?
- Do we need to ditch spurious modes before calculating pseudospectrum? That tearing instability paper did this in a 
  neat way: basically projected out the spurious subspaces from the operator before calculating pseudospectrum
- Consider presenting at APS DPP?
  - Abstract deadline July 15
  - Possibly present alongside weird modes that appear at strong magnetic fields + low Rm 
    (I've found this in sinusoidal shear flow, haven't looked in tanh layer yet)