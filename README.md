# nonmodal-MHD-KH
Exploring nonmodal effects in KH-unstable shear flows in MHD using 
[eigentools](https://github.com/DedalusProject/eigentools), a set of tools for studying eigenvalue problems that's 
built on top of [Dedalus](dedalus-project.org). For two great examples of why one should care about these issues, see 
[Kaminski et al. 2014](https://doi.org/10.1017/jfm.2014.552) and 
[Kaminski et al. 2017](https://doi.org/10.1017/jfm.2017.396), which consider KH-unstable shear layers in the 
hydrodynamic case with density stratification.


This is an on-going, collaborative project. We're sharing our progress publicly for now in the spirit of open science 
and communication of results, but we ask that you please reach out to talk about attribution/acknowledgements/citation 
if you use any code, ideas, or results you find here in your own work.

The equations we're working with are derived in `tex/equations.tex`. A script for calculating the epsilon-pseudospectra 
of that system is given in `python/MHD_OrrSomm_pseudospec.py`. `python/MHD_OrrSomm_EVP.py` is a starting point for 
solving the eigenvalue problem without all the great tools provided by eigentools, and can probably be ignored. 
Various internal notes are in `internal_notes_logs/`.