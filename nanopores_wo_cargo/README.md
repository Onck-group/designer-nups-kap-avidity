

# Nanopore simulations - without cargo

The general runscript is found in `jobscript` , including all GROMACS commands. Simulations were carried out with GROMACS 2018.

Relevant files to re-run the GROMACS simulations, such as forcefield files, tabulated potentials, protein structures, etc. are in the respective sub-folders, labeled per NupY-variant.

Note that various equilibration steps are necessary to relax the structure of the nanopore assemblies. The resultings structures (`hyper_1.pdb`, `hyper_2.pdb`, etc.) are included to skip these steps in reproduction.

