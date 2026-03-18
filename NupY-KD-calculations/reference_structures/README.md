This respository contains an old method of generating the `.pdb`-files for the single NupY-proteins (extracted from a nanopore environment). It is not used actively / necessary for reproducing the results of this study.



Some notes:
-A preceding GP-peptide is added to the sequence, and an additional Cys-residue is added to the end
-Three constructs differ very slightly from the sequences used in the experiments. These are C/H=0.1912+0.3824, dFG=13 (one, resp. 2 additional FG-motifs in the sequences used in the simulations, corrected in the genetic construct for the experiments), and C/H=0.5736, dFG=7 (one T and K residue have been swapped due to issues with dyad symmetries).
-C/H=0.0239 has not been studied using selectivity simulations, and none of the C/H=0.0239 constructs is studied experimentally. This protein is incredibly cohesive, and we do not expect it to lead to any feasible results if we try to purify it. 
-C/H=0.0239, dFG=7 is not included since none of the 10^4 preliminary designs displayed sufficient disorder.
- `.pdb`-files were extracted from a nanopore simulation. These used relatively old codes, which do not properly include the residue numbers and residue names. Consider generating a fresh topology + structure instead (and perhaps replace coordinates accordingly).

