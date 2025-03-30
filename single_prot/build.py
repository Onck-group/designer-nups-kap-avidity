###################################################################
# 
# 
#
# Originally written by A. Jansen, MSc, edited by H.W. de Vries (2020)
#
# Uses sys, np, os python packages, and any installed version of
# GROMACS (i.e. gmx) to pre-process structure files.
#
# Creates all necessary files to run a coarse-grained MD/SD simulation of
# a disordered protein using the 1-BPA forcefield. 
# See Ghavami et al. BPJ 2015
# See Ghavami et al. JCTC 2013.
#
# Input: 
# Protein name (str) = the name of the protein. The sequence of the protein
# should be present in the same directory, in a file that has the same name.
# The file itself should be formatted such that the first line comprises
# '%ProteinName', and the second line contains the full sequence (so: no FASTA!)
# 
# Output: 
# Creates a directory <protein name> and copies all files from the folder
# grom to this directory. Also, the following files are created:
#
# chain_org.pdb + chain_org2.pdb = structure files with the 1-BPA resolution
# 								   protein chain, fully extended. The second
#								   version of the .pdb file contains box vectors.
# proteinname.itp = topology file (GROMACS) for the individual protein.
# Forcefield.itp = forcefield file with 1-BPA parameters.
# chain.ndx = bookkeeping file with index numbers and several pre-written groups.
#             This is used for analysis, and can be edited post-simulation.
# chain.top = gromacs-readable topology file, which loads all other topology files
#             and indicates the occurrences of all system components
#
#
# Note that customized tables (with the LJ and coulombic potentials), the custom
# bonded potentials, etc. are all copied from the grom-folder. In case a new table.xvg
# file needs to be made (e.g. the ionic strength changed), then re-make the table file
# using generatetable.py
###################################################################
import sys
import numpy as np


name      = str(sys.argv[1])

# PARAMETERS for the 1-bpa forcefield - shouldn't be touched.

BBexc     = 3		# 				Number of neighbour exclusions.
bl        = 0.38	#nm 			Bond length between two beads (0.38 is carbon-carbon bond length).
K_b	      = 8038  	#kJ/nm^2/mol 	Spring constant of said bond.

nbfunc	  = 1		# 1 = Use LJ-like potential (c6, c12). 2 = Use Buckingham-like potential. Set this to 1.0.
comb_rule = 1		# Mixing rule. 1 = intepret c6 and c12 directly. Set this to 1.0.
gen_pairs = 'no'	# Should be set to no.
fudgeLJ   = 1.0		# Hoax factor. Set this to 1.0.
fudgeQQ   = 1.0		# Hoax factor. Set this to 1.0.

eps_HP    = 13			#kJ/mol
eps_rep   = 10			#kJ/mol
alpha     = 0.27		#Free parameter

mass      = 124			#Da 			Average residue mass.
sigma     = 0.60		#nm 			Average bead radius.
catpi_eps = sys.argv[2]
charge    = {'A':0,    'R':1,     'N':0,    'D':-1,    'C':0,    'Q':0,    'E':-1,    'G':0,    'H':0,    'I':0,    'L':0,    'K':1,     'M':0,    'F':0,    'P':0,    'S':0,    'T':0,    'W':0,    'Y':0,    'V':0}
HP_dict   = {'A':0.70, 'R':0.005,     'N':0.33, 'D':0.005,'C':0.68, 'Q':0.675, 'E':0.005,'G':0.41, 'H':0.53, 'I':0.98, 'L':1.00, 'K':0.005,'M':0.78, 'F':1.00, 'P':0.65, 'S':0.45, 'T':0.51, 'W':0.96, 'Y':0.82, 'V':0.94}

###################################################################
# MODULES/GENERAL FUNCTIONS
###################################################################

import os
import numpy as np

def getstring(file,line,column):
	'''
	Assumes that the input file contains:
	 on line 1 %<name> 
	 on line 2 <sequence>
	'''
	x = open(file,'r')
	for y, z in enumerate(x):
		if y == line-1:
			string = z.split()[column-1]
	x.close()
	return string



# Dictionaries for translating to full residue names
aa_lib = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

res_lib = {'C' : 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
           'I' : 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN',
           'G' : 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP',
           'A' : 'ALA', 'V': 'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}


###################################################################
# MAKE DIRECTORY STRUCTURE
###################################################################

# WORKING DIR

os.system('echo Working dir is ${PWD}')

# CLEANUP


if os.path.isdir('Output/%s' % (name)) is True:
	print('Removing previous output dir...')
	os.system('rm -r Output/%s' % (name))

# LOAD PROTEIN SEQUENCE

print('Reading protein sequence...')

seq = getstring(name, 2, 1)
res = [res_lib[aa] for aa in seq]
print(res)

# CREATE NEW OUTPUT DIR

print('Creating output directory...')

os.makedirs('Output/%s' % (name))

# COPY GROM FILES

print('Copying grom files...')

os.system('cp -r ./grom/* ./Output/%s/' % (name))

###################################################################
# CREATE chain.top
###################################################################

print("Building chain.top...")

f = open('chain.top','w+')

f.write('[ defaults ]\n')
f.write('; nbfunc  comb-rule  gen-pairs	 fudgeLJ  fudgeQQ\n')
f.write('    %s        %s           %s	   %s      %s\n' % (nbfunc, comb_rule, gen_pairs, fudgeLJ, fudgeQQ))
f.write('\n')
f.write('#include \"forcefield.itp\"\n')
f.write('#include \"%s.itp\"\n' % (name))
f.write('\n')
f.write('[ system ]\n')
f.write('; Name\n')
f.write('%s\n' % (name))
f.write('\n')
f.write('[ molecules ]\n')
f.write('; Compound	#mols\n')
f.write('%s		1\n' % (name))

f.close()

os.system('mv chain.top ./Output/%s' % (name))

###################################################################
# CREATE .itp
###################################################################

print("Building %s.itp..." % (name))

# REBUILD CHAIN

angles = []; dihedrals = []; rebuild = []

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

# WRITE HEADER

f = open('%s.itp' % (name),'w+')

f.write('[ moleculetype ]\n')
f.write('; Name      nrexcl\n')
f.write('%s	    %s\n' % (name, BBexc))
f.write('\n')

# WRITE [ ATOMS ]

f.write('[ atoms ]\n')
f.write('; nr type resnr residue atom cgnr\n')

for i in xrange(0,len(seq)):
	f.write('%s \t %s \t %s \t %s \t %s \t %s\n' % (i+1, seq[i], i+1, res[i], seq[i], i+1))

# WRITE [ BONDS ]

f.write('\n[ bonds ]\n')
f.write('; i  j funct   c0       c1\n')

for i in xrange(1,len(seq)):
	f.write('%s \t %s \t %s \t %.3f \t %.3f\n' % (i, i+1, 1, bl, K_b))

# WRITE [ ANGLES ]

f.write('\n[ angles ]\n')
f.write('; ai    aj    ak     funct  table_number  k(kJ/mol)\n')

for i in xrange(3,len(seq)+1):
	f.write('%s \t %s \t %s \t %s \t %s \t %s \n' % (i-2, i-1, i, 8, angles[i-3], 1))

# WRITE [ DIHEDRALS ]

f.write('\n[ dihedrals ]\n')
f.write('; ai    aj    ak    al    funct   table_number   k(kJ/mol)\n')

for i in xrange(4,len(seq)+1):
	f.write('%s \t %s \t %s \t %s \t %s \t %s \t %s\n' % (i-3, i-2, i-1, i, 8, dihedrals[i-4], 1))

f.write('\n')
f.close()

os.system('mv %s.itp ./Output/%s' % (name, name))

###################################################################
# CREATE chain_org.pdb
###################################################################

print('Building chain_org.pdb...')

# CREATE PROTEIN COORDINATES

x_cycle = [0.000, 3.800, 5.700, 3.800, 0.000, -1.9000]
y_cycle = [0.000, 0.000, 3.291, 6.582, 6.582, 3.291]

x_coords = []
y_coords = []
z_coords = []

for i in xrange(0,len(seq)):
	x_coords.append(x_cycle[i % 6])
	y_coords.append(y_cycle[i % 6])

val = 0
for i in xrange(1,len(seq)+2):
	z_coords.append(val)

	if i % 6 == 0:
		val += 5

z_coords.remove(5)

# WRITE HEADER

f = open('chain_org.pdb','w+')

f.write('TITLE     CHAIN\n')
f.write('MODEL        1\n')

# WRITE COORDINATES

for i in xrange(0,len(seq)):
	f.write('{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n'.format('ATOM', i+1, seq[i], '', res[i], '', i+1, '', x_coords[i], y_coords[i], z_coords[i], 1.00, 0.00))

# WRITE CONECTS

for i in xrange(1,len(seq)):
	f.write('{:6s}{:5d}{:5d}\n'.format('CONECT', i, i+1))

# f.write('{:6s}{:5d}{:5d}\n'.format('CONECT', len(seq), 2*len(seq)))

f.write('TER\n')
f.write('ENDMDL')

f.close()

# CHANGE TO PERIODIC BOX

box_size = np.ceil(len(seq)/12.0) + 12

os.system('gmx editconf -f chain_org.pdb -o chain_org2.pdb -box %s %s %s -center %s %s %s &> /dev/null' % (box_size, box_size, box_size, box_size/2.0, box_size/2.0, box_size/2.0))
os.system('mv chain_org.pdb ./Output/%s' % (name))
os.system('mv chain_org2.pdb ./Output/%s' % (name))

# ###################################################################
# CREATE forcefield.itp
# ###################################################################

print('Building forcefield.itp...')

f = open('forcefield.itp','w+')

# WRITE HEADER

f.write('[ atomtypes ]\n')
f.write(';name   mass  charge ptype        c6             c12\n')

# WRITE SAME-SAME ENERGIES

A1 = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

for i in A1:
	c6  = eps_HP * ( (HP_dict[i]*HP_dict[i])**alpha )**0.5
	c12 = eps_rep

	f.write('%s \t %s \t %s \t %s \t %s  \t %E  \t %E \n' % (i, i, mass, charge[i], 'A', c6, c12))
f.write('\n[ nonbond_params ]\n')
f.write('; i     j   func        C            A\n')

for i in xrange(0, len(A1)):
        for j in xrange(i, len(A1)):
                if A1[i] != A1[j]:

                        # CAT-PI
                        if A1[i] == 'R' and A1[j] in ['F','Y','W']:
                                c6  = float(catpi_eps)
                                c12 = float(catpi_eps)

                                f.write('%s \t %s \t %s \t %E \t %E \n' % (A1[i], A1[j], 1, c6, c12))
                        if A1[i] == 'K' and A1[j] in ['F','Y','W']:
                                c6  = 0.33* float(catpi_eps)
                                c12 = 0.33* float(catpi_eps)

                                f.write('%s \t %s \t %s \t %E \t %E \n' % (A1[i], A1[j], 1, c6, c12))                                

f.write('\n')
f.write('\n')
f.close()

os.system('mv forcefield.itp ./Output/%s' % (name))

# ###################################################################
# CREATE chain.ndx
# ###################################################################

print('Building chain.ndx...')

f = open('chain.ndx','w+')

f.write('[SYSTEM]\n')

count = 0
for i in xrange(0, len(seq)):
	f.write('%s  ' % (i+1))

	count += 1
	if count == 300:
		f.write('\n')
		count = 0

f.write('\n')


f.write('\n\n[AA]\n')

count = 0
for i in xrange(0, len(seq)):
        if seq[i] not in ['R','K','F','Y','W']:
                f.write('%s  ' % (i+1))

                count += 1
                if count == 300:
                        f.write('\n')
                        count = 0

f.write('\n\n[RK]\n')

count = 0
for i in xrange(0, len(seq)):
        if seq[i] == 'R' or seq[i] == 'K' :
                f.write('%s  ' % (i+1))

                count += 1
                if count == 300:
                        f.write('\n')
                        count = 0

f.write('\n\n[FYW]\n')

count = 0
for i in xrange(0, len(seq)):
        if seq[i] in ['F','Y','W']:
                f.write('%s  ' % (i+1))

                count += 1
                if count == 300:
                        f.write('\n')
                        count = 0

f.write('\n')



f.close()

os.system('mv chain.ndx ./Output/%s' % (name))

######################################################################################################################################
# CREATE table.xvg
######################################################################################################################################

###################################################################
# SHARED FUNCTIONS
###################################################################

