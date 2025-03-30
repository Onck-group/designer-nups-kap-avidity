#!/bin/bash
gmx grompp -f hyper_EM.mdp -p hyper.top -o hyper_EM2.tpr -c hyper_EQ.pdb -n sys.ndx -maxwarn 4
srun gmx_mpi mdrun -s hyper_EM2.tpr -o hyper_EM2.trr -c hyper_EM2.pdb -g EM2.log -e em2.edr -v -tableb table_a*.xvg table_d*.xvg

gmx grompp -f hyper_MD.mdp -p hyper.top -o hyper_MD.tpr -c hyper_EM2.pdb -n sys.ndx -maxwarn 4
srun gmx_mpi mdrun -s hyper_MD.tpr -o hyper_MD.trr -c hyper_MD_out.pdb -g MD.log -e md.edr -v -tableb table_a*.xvg table_d*.xvg  -dlb yes


