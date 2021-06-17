# June 17th, 2021
### Progress
- When pushing to higher Re, the subspace dimension ("k" in the documentation) turns out to be super important
- Lots of striking pseudospectra, though

### Over the next 2 weeks
- Write script to save spectra/pseudospectra/eigenvectors to .hdf5 files, rather than just saving plots
- Go back and play with "k" more carefully, including in Orr-Sommerfeld problem
- Clean up sinusoidal shear flow problem from DDC work, leverage the fact that it's a non-generalized eigenvalue 
problem, meaning we can calculate the pseudospectra the old-fashioned, dense way
- Stop plotting low and high eigenvalues! Just plot good eigenvalues

# June 3rd, 2021
### Progress
- Calculated my first pseudospectrum! The functioning script is in nonmodal-MHD-KH/python/MHD_OrrSomm_pseudospec.py
  
### Questions/issues/discussion
(Strikethrough things we resolved)
- ~~Dissipative or dissipationless spectra?~~
- Do you usually calculate pseudospectra at all unstable wavenumbers? Most-unstable wavenumber?
- ~~Do we need to ditch spurious modes before calculating pseudospectrum? That tearing instability paper did this in a 
  neat way: basically projected out the spurious subspaces from the operator before calculating pseudospectrum~~
- ~~What's next?~~
  - ~~Come up with a relevant metric characterizing pseudospectra, then~~ do parameter scans
  - ~~3D perturbations?~~
  - ~~Linear optimal perturbations?~~
- Can we talk about how Eigentools calculates pseudospectra under the hood?
  - ~~Noticed calc_ps calls EP.solve(sparse=True). If I first call EP.solve(sparse=False), check the % difference between 
    evalues_low and _high for the unstable mode, then check it again after calling calc_ps, the % difference goes up,
    causing the mode to go from not-rejected to rejected.~~
- ~~Consider presenting at APS DPP?~~
  - ~~Abstract deadline July 15~~
  - ~~Possibly present alongside weird modes that appear at strong magnetic fields + low Rm 
    (I've found this in sinusoidal shear flow, haven't looked in tanh layer yet)~~
- I usually add figure files (*.pdf, etc) to .gitignore. What's the strategy when working with Overleaf+github?

### Post-meeting recap
- Crank up the Reynolds number. In stratified case, things get more interesting at Re ~ 10^3 - 10^4 and non-modal 
  growth increases
- More parameter explorations before moving on to 3D or optimal perturbations
- Jeff looking into MacTaggart paper, projecting out spurious subspaces idea
- Nothing wrong with dissipationless case
- APS DPP sounds fine, provided understanding that it's ok to present works in progress there
- Talked about available computing resources
- Minor differences between sparse/dense solves expected