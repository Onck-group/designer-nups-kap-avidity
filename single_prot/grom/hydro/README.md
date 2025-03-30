## Running HYDRO++ on IDP simulations

HYDRO++ operates on any input .pdb file, and is applicable to systems comprising only Ca atoms. Run parameters were adapted from Ghavami et al. JCTC 2013 and Ghavami et al. BPJ 2014. 

To perform the analysis, take the following steps:

- Convert a subset of trajectory frames to the input .pdb file, after removing pbc using: 
  - `echo 0 | gmx trjconv -f MD_chain.trr -s MD_chain.tpr -pbc nojump -o hydro.xtc`
  - `echo 0 | gmx traj -s MD_chain.tpr -n chain.ndx -f hydro.xtc -oxt MD_chain_hydro.pdb -b 500000`. The analysis in this case starts from 500 ns on (and uses the subsequent 1 us of data).
- Call `./runhydro.py` via the command line (make sure that you copied the `runhydro.py` script to the same file as the `MD_chain_hydro.pdb` file!)
- Output is written in `hydro.out`

Note that HYDRO++ is written in FORTRAN and is able to run in parallel on computer clusters as well (this might speed things up notably).

