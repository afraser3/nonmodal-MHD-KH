# nonmodal-MHD-KH
Exploring nonmodal effects in KH-unstable shear flows in MHD using 
[eigentools](https://github.com/DedalusProject/eigentools), a set of tools for studying eigenvalue problems that's 
built on top of [Dedalus](dedalus-project.org). 
This work uses v2 of Dedalus, not the latest "d3". Installation instructions for Dedalus and eigentools will not be 
provided here, but are provided in the documentation for those packages.

This work has been accepted for publication in Physical Review E at [this DOI](https://doi.org/10.1103/vsm7-1zj4). 
A preprint can be found [here](https://arxiv.org/abs/2510.03141). The purpose of this repository is to share the 
necessary source files for generating the data and the figures in that publication.

A lot of the files in this repo are outdated / "legacy" in the sense that they are not used in generating the figures
that went into the final publication. They are retained in this repo just in case they are somehow useful to others
in the future. 

The overall workflow for generating data that went into the publication is as follows:
- Run `python3 scan_linear_optimals_3d.py PATH/TO/CFG_FILE.cfg` to calculate and save the LOPs for a given set of 
parameters across a (specified) range of wavenumbers; see `run_scripts/pleiades_LOPs_3D.pbs` for an example that runs
this on a cluster using PBS.
- For the list of specific CFG files that went into our publication, inspect `run_dict` and `run_dict_pleiades` in 
`paper_figure_plotters/helper_functions.py`. These CFG files can be found in `LOPs_3D/config_files/`
- Once you've generated LOPs for each of these CFG files, run `paper_figure_plotters/Gs_kopt_LOPs_IVPs_smallbox.py` and 
`paper_figure_plotters/Gs_kopt_LOPs_IVPs_bigbox.py` to create the figures we used in the paper

You should now have the figures in the paper.