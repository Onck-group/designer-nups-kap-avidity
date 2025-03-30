## Gene optimization of NupY variants

Gene optimization was typically carried out using the following steps:

1. Translation of the designed 800-residue NupY mutant to an E.Coli (K12 strain) gene using the NovoProLabs tool. No restriction enzymes were excluded originally, but consider doing this at this stage (or check manually later). Tool is available at: https://www.novoprolabs.com/tools/codon-optimization
2. Addition of the His8 tag, and a 3C protease cleave site (N-end) and the addition of a cysteine residue for anchoring (C-end) together with the respective codons and a stop-codon. For this, we use the python-code `fasta_converter.py`. 
    - For this we typically use ATG GGC CAC CAT CAC CAT CAC CAT CAC CAT GAT TAC GAT ATT CCA ACG ACC CTG GAA GTT CTG TTC CAG GGG CCC , which codes for M   G   H   H   H   H   H   H   H   H   D   Y   D   I   P   T   T   L   E   V   L   F   Q   G   P
3. At this stage, we optimize the CAI for expression in E.Coli and removal of rare codons. For this we use the tool http://gcua.schoedl.de/sequential_v2.html with E.Coli K12 strain. The easiest way to proceed is to let the output be e-mailed as a .txt file, and afterwards run piped bash commands: a grep-search for the words 'grey' and/or 'red' to find the rare codons, and then use 'awk' to print the residue numbers. These can be pasted into any text editor, where new-line symbols can be removed by ',<space>' to form a comma-separated list of residues that need to have their codons swapped out for less-rare options. 
4. For the above we use an in-house code (`optimize_gene_pointwise.py`) that is able to take amino acid numbers or base numbers to perform codon swaps. It uses an adapted codon table that is biased towards common residues in E.Coli. The AA numbers can be included in the `aa_numbers_to_change` list, after which a simple list comprehension converts these to the relevant codon numbers in the `optimize_gene_pointwise.py` code.
5. After the CAI optimization we perform an analysis of the dyad symmetry (i.e. highly energetic stem loops/hairpins) and check for presence of rho-independent terminators. 
    - For this we use the tools ARNOLD and FindTerm (latter can only be used 20 times/day per IP address). 
    - ARNOLD requires a FASTA sequence can be found at : http://rssf.i2bc.paris-saclay.fr/toolbox/arnold/ , where we used the standard settings (scan in both directions). 
    - For FindTerm, the box for display of all putative terminators should be ticked, and the energy threshold is set at -1.5 units. http://www.softberry.com/berry.phtml?topic=findterm&group=programs&subgroup=gfindb&advanced=on
    - The above produces a range of residues in which terminators exist. This range is fed to the `optimize_gene_pointwise.py`-script, which performs CAI-optimizing mutations, after which we repeat the checks until satisfactory.
6. Further checks for stem loops/hairpins can be made using CloneManager (proprietary software). Accidental restriction enzyme sites can be removed using the same mutation script. Note that in some cases the use of optimized codons is not possible, due to the limited number of alternative codons (for example, when a terminator or highly energetic stem loop persists). Thee `alternative=False` line can be set to `True` in `optimize_gene_pointwise.py` to use alternative codons without affecting the overall CAI too much.
7. The whole is cloned into pET28-a with NcoI and Xhol restriction enzyme sites.

Note that for NupY mutants, we re-used the originally optimized sequence for the extended domain. For this, the `replace_ext_domain.py` script is available.

The resulting sequence (AA) is: GP-designed sequence-C, totaling 803 residues.

### Known issues:
- The original output of `convert_fasta.py` cannot be easily read into the `optimize_gene_pointwise.py` script (the latter assumes a single-line formatted file). Some additional commands are necessary, which can be found in an alternative version: `optimize_gene_pointwise_alt.py`
- Restriction enzymes should be excluded in the novopro-design step as well, ideally. Alternatively, ensure that these are not present on erroneous positions.
- Application of the CAI optimization step shouldn't involve swapping of codons in the range of codons 1:75, this is specifically pre-optimized already
- Re-use of the original NupY genetic sequence for the extended domain (residues 611-800. or 613-801) might lead to a feature appearing in the findTerm programme's output (near 1951). However, this part of the protein was adequately expressed in earlier experimental work and is not expected to alter expression.
- Use the `alternative=True`-function cautiously, if it is possible to get rid of stem loops and terminators without it, then this is preferred.

