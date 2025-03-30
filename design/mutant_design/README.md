### Design of NupY mutants
The codes in this folder were used to create variants of NupY, where the cohesive domain (residues 1-610) were varied in a controlled but systematic way, such that properties like the FG-spacing or C/H ratio were varied.


- `make_FG_spacing_mutant.py` : this code considers a basefile (NupY50.txt) and works in various steps. First, the code swaps out all the present FG/GLFG motifs for other motifs that have the same hydrophobicity on average, according to the epsilon-i values used in the 1-BPA model (see Ghavami et al. BPJ 2014). These pairs are selected using `test_comobos.py`. Notably, aromatic residues and combinations of two order-promoting residues are avoided. Following the replacement of the FG-motifs, FG and GLFG motifs are inserted on new positions, based on the chosen spacing. Any spacer residues that are already present on sequence indices where FG/GLFG motifs are now inserted, are stored. Then, these former spacer residues replace residues that are in FG/GLFG-replacing motifs (i.e. the original FG/GLFG positions). This preserves the C/H ratio and minimizes the alterations to the sequence.
- `make_CH_mutant.py` : this code considers a basefile (i.e. NupY or a FG-spacing mutant) and mutates spacer residues (non FG/GLFG-motifs) such that the C/H ratio of the protein converges towards the target value. The signed net charge per residue is enforced to be within 2.5% of the value in the cohesive domain of the ..


- `determine_disorer_scores.m` processes the output of SPOT-Disorder-Single, and ranks designs with fully disordered cohesive domains based on their average disorder score. This can be used to make a selection from the pool of designs.
