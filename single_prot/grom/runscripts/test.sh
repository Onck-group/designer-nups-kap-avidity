#!/bin/bash

gmx grompp -f EM_chain.mdp -c chain_org2.pdb -p chain.top -n chain.ndx -o EM_chain.tpr -maxwarn 4
gmx mdrun -s EM_chain.tpr -o EM_chain.trr -c EM_chain.pdb -g EM_chain.log -e EM_chain.edr -v -tableb table_a*.xvg table_d*.xvg -maxh 0.5

gmx grompp -f EQ_chain.mdp -p chain.top -o EQ_chain.tpr -c EM_chain.pdb -n chain.ndx -maxwarn 4
gmx mdrun -s EQ_chain.tpr -o EQ_chain.trr -c EQ_chain.pdb -g EQ_chain.log -e EQ_chain.edr -v -tableb table_a*.xvg table_d*.xvg -maxh 0.5

gmx grompp -f MD_chain.mdp -p chain.top -o MD_chain.tpr -c EQ_chain.pdb -n chain.ndx -maxwarn 4
gmx mdrun -s MD_chain.tpr -o MD_chain.trr -c MD_chain.pdb -g MD_chain.log -e MD_chain.edr -v -tableb table_a*.xvg table_d*.xvg 


