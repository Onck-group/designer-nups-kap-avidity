import MDAnalysis as mda 
import sys

import numpy as np
import os

#THIS SCRIPT ASSUMES THAT YOU HAVE A VERSION OF GROMACS INSTALLED OR LOADED.


for i in range(5):
	print(sys.argv[i])


proteinFile=sys.argv[1]
proteinCode,extension = proteinFile.split(".",1)

#The protein also needs a name for further usage. Insert it as the second argument.
proteinName=sys.argv[2]
#prepare a directory for output.
outputDir = './'

if sys.argv[3] == 'True':
	top_switch = True
else:
	top_switch = False

if sys.argv[4] == 'True':
	pdb_switch = True
else:
	pdb_switch = False



#Three-letter to one-letter conversion dictionary (BioPython could also do this but I'm lazy)
aa_lib = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

res_lib = {'C' : 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
           'I' : 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN',
           'G' : 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP',
	       'A' : 'ALA', 'V': 'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET',
	       'F1' : 'PHE', 'G1' : 'GLY'}





def prepare_protein(proteinFile, proteinCode, proteinName ):

	#Prepare the protein: 
	#1. take out all the lines that contain a C-alpha atom ("CA") s
	#2. write it to a separate file.
	protein_object=mda.Universe(proteinFile)
	protein_coordinates = protein_object.select_atoms('all')
	sequence_length = len(protein_coordinates)
	aa_sequence = []
	for i in range(sequence_length):
		cur_res_3 = protein_coordinates[i].resname
		cur_res_1 = protein_coordinates[i].name 
		if i < sequence_length-1:
			next_res  = protein_coordinates[i+1].name
			prev_res = protein_coordinates[i-1].name
		aa_sequence.append(cur_res_1[0])

		if cur_res_1 == 'F' and next_res == 'G':
			protein_coordinates[i].name = 'F1'
		if cur_res_1 == 'G' and prev_res == 'F1':
			protein_coordinates[i].name = 'G1'


	return protein_object, protein_coordinates, aa_sequence


def rebuild_protein_topology(aa_sequence, proteinName, protein_coordinates):
# This function writes a topology file for an input folded protein.
# It assigns the relevant backbone potentials (based on whether the AA's are G (flexible), P (stiff) or X (anything else in between)).
# The first part of this function is mainly based on earlier work from A. Ghavami and the python adaption by A. Jansen. 
# the function also generates a network of harmonic bonds to preserve secondary and tertiary structure. 
# The bonded network employs two different cut-offs.
	seq=aa_sequence[:]
	res=[res_lib[aa] for aa in seq]
	angles = []; dihedrals = []; rebuild = []
	BBexc=3
	bl=0.38
	K_b=4019.328*2
	for i in seq:
		if i == 'G':
			rebuild.append(1)
		elif i == 'P':
			rebuild.append(2)
		else:
			rebuild.append(3)
		for i in xrange(1,len(rebuild)-1):
			if rebuild[i+1] != 2:
				angles.append(rebuild[i])
			else:
				angles.append(int('%s%s' % (rebuild[i],rebuild[i+1])))

		for i in xrange(1,len(rebuild)-1):
			dihedrals.append(int('%s%s' % (rebuild[i],rebuild[i+1])))
	
	
	f = open('%s.itp' % (proteinName),'w+')

	f.write('[ moleculetype ]\n')
	f.write('; Name      nrexcl\n')
	f.write('%s	    %s\n' % (proteinName, BBexc))
	f.write('\n')
	
	f.write('[ atoms ]\n')	#Begin with the definition of all the atoms
	f.write('; nr type resnr residue atom cgnr\n')
	for i in xrange(0,len(seq)): #TK: REPLACE THIS WITH PROTEIN NAME indices!
		f.write('%s \t %s \t %s \t %s \t %s \t %s\n' % (i+1, protein_coordinates[i].name, i+1, res[i], protein_coordinates[i].name, i+1))

	f.write('\n[ bonds ]\n') #First topology section: bonds
	f.write('; i  j funct   c0       c1\n')
	for i in xrange(1,len(seq)):
		f.write('%s \t %s \t %s \t %.3f \t %.3f\n' % (i, i+1, 1, bl, K_b))

	f.write('\n[ angles ]\n') #Second section: angles
	f.write('; ai    aj    ak     funct  table_number  k(kJ/mol)\n')
	for i in xrange(3,len(seq)+1):
		f.write('%s \t %s \t %s \t %s \t %s \t %s \n' % (i-2, i-1, i, 8, angles[i-3], 1))

	f.write('\n[ dihedrals ]\n') #Third section; dihedral potentials
	f.write('; ai    aj    ak    al    funct   table_number   k(kJ/mol)\n')
	for i in xrange(4,len(seq)+1):
		f.write('%s \t %s \t %s \t %s \t %s \t %s \t %s\n' % (i-3, i-2, i-1, i, 8, dihedrals[i-4], 1))
	f.write('\n')
	
	f.close		






protein_object, protein_coordinates, aa_sequence = prepare_protein(proteinFile, proteinCode, proteinName)

if top_switch == True:
	rebuild_protein_topology(aa_sequence, proteinName, protein_coordinates)

if pdb_switch == True:

	protein_coordinates.write('%s_rebuilt.pdb' % (proteinName))




