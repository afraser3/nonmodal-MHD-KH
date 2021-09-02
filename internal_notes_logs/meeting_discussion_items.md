# September 2nd, 2021
### Progress
- Worked up to Re = 2500
- Messed around with "mu" parameter to get better-behaved pseudospectra at high MA
- Took time off! And started fellowship proposals.

### Recap
- Everyone is on-board with including this stuff in fellowship proposals, with the understanding that we intend to keep collaborating either way because this project is great
- Noticed that it's not just weird behavior in pseudospectra, but the spectrum might be funny, mainly weird omega = 0, gamma < 0 modes e.g. in run AF. Let's plot those modes and see what's going on.
- Created Google Drive folder to dump run files (.h5, .png, and .pdf) into and share with team
- Should definitely create overleaf document to maintain organization of ideas/results as this stretches on
- Immediate next steps, maybe for DPP:
  - To help diagnose this "mu" and "k" business in weird pseudospectra behavior, take the linear operator for Kolmogorov flow (which I've already coded up as an ordinary EVP in 2D, 3D, and for streaky flows, due to ongoing work w/ P. Garaud) and calculate the pseudospectra directly using resolvent, then compare with pseudospectra obtained with this new method
    - This is too technical and uninteresting to deal with before DPP, but we're gonna do this before publishing
  - Pseudospectra of 3D modes
    - Easy! Do this before DPP if the following doesn't pan out
  - Linear optimals
    - Conceptually manageable, but it's anybody's guess as to how long it will take to converge on a solution for a single parameter point, thus who knows how long parameter scans will take
    - Following Kaminski PhD Thesis Ch 2, solve for the adjoint equations, quickly code them up in Dedalus, and just see how long a single parameter point takes
    - For this test, we care more about how long it takes than the actual result, so don't bother with careful selection of norm
    - If it takes long, kick this down the road for post-DPP. If it's doable, let's include this at DPP.
    - Sure would be fun to do nonlinear optimals, compare result of energy norm vs dissipation norm, and see how much stable mode excitation is in initial condition or occurs along the way to the optimized time.

# July 1st, 2021
### Progress
- Separate scripts for calculating and saving pseudospectra vs plotting from file
- Not much news on figuring out important "k". Importance of larger k almost seems 
  to track with importance of larger Nz? As in, at higher Re.
  
- Insufficient k always makes pseudospectra seem smaller/more trivial

### Questions/issues
- How to do calc_ps with MPI?
- What should we aim to finish in time for DPP abstract submission?
- Things it would be nice to say:
  - Which regions of parameter space, broadly, have more non-normality?
  - When KH is barely stabilized by B, is there significant potential for non-modal growth?
  - Maybe: do the different varieties of MHD KH have significantly different degrees of non-normality?
  
- What can we reasonably promise to finish by Nov?
  - Linear optimal perturbations, comparison of transient vs modal growth?
  - Maybe nonlinear simulations of linear optimals?

### Recap
- DPP:
  - Nonlinear DNS and/or optimals would be great to do, but too ambitious for this year. 
    Linear optimals + pseudo is plenty.
  - Don't sweat coming up with new results in the next week just for sake of abstract. We have enough for an abstract.
  - Will send around abstract draft next week. Note pseudospectra/non-normality more familiar at DFD than DPP, use jargon accordingly.
  
- Pseudospectra only getting more interesting as k increases is a relief
- Jeff pointed to [this link](https://github.com/jsoishi/mri_prefers_3d/blob/715022a2f8ba8361bac7cdcaeb52e70ab2257468/python/mri_single_yz_mode.py#L81)
  and, even better, [this link](https://github.com/jsoishi/mri_prefers_3d/blob/715022a2f8ba8361bac7cdcaeb52e70ab2257468/python/mri.py#L230),
  for examples of MPI plus eigentools
  
- Don't sweat doing a kx scan at every parameter point. Maybe look into it here and there, but fixing kx and scanning other stuff is probably fine
- Now that saving to hdf5 is set up, once MPI is set up, can just scan lots of MA/Re/Pm at once and at sufficiently high k/Nz on Lux (local cluster)

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