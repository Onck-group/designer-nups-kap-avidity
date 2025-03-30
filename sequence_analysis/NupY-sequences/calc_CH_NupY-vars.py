import sys
import numpy as np
import random
from math import ceil
import copy
import os


def getseq(file,line,column):
	'''
	reads in the sequence from a source file, assuming that everything is placed on one very long line.
	This might also work for FASTA etc, but requires a different way of being called. Use with caution!
	Credits to Anton Jansen for putting this black box together.
	'''
	x = open(file,'r')
	for y, z in enumerate(x):
		if y == line-1:
			string = z.split()[column-1]
	x.close()
	list_string=[aa for aa in string]
	return list_string

def pick_scale(scale_name='ghavami'):
	'''
	pick_scale(scale_name)
	returns a hydrophobicity scale, depending on the chosen input.
	Options are 'ghavami', 'ww', 'kd', 'eis', which correspond to the following:
		'ghavami' = scale used in the 1-BPA MD model
		'ww'      = Whimley-White hydrophobicity scale
		'kd'      = kyte-doolittle hydrophobicity scale
		'eis'     = eisenberg hydrophobicity scale. Eisenberg 1984
		'count'   = count a hydrophobic residue (AIFLWYV)
	note that the output has to be CONSISTENT with the order in which
	the different AA's are presented in the original containing vector.
	The values, normalized between 0 and 1, were obtained by: (scale+|min(scale)|)/max(scale).
	'''
					 #Ala,Arg,Asn,Asp,Cys,Gln,Glu,Gly,His,Ile,Leu,Lys,Met,Phe,Pro,Ser,Thr,Trp,Tyr,Val
	amino_acid_names=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
	#    1-based       1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20
	#    0-based       0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10, 11, 12, 13, 14, 15, 16, 17, 18, 19

	if scale_name == "ghavami":
								#   ['A','R','N' , 'D'  , 'C', 'Q', 'E'  , 'G', 'H','I','L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
		
		hydrophobicities = np.array([0.7,0.005,0.33,0.005,0.68,0.64,0.005,0.41,0.53,0.98,1,0.005,0.78,1.00,0.65,0.45,0.51,0.96,0.82,0.94])
		return amino_acid_names, hydrophobicities

	elif scale_name == "ww":
		hydrophobicities = np.array([0.48,0.31,0.41,0.20,0.58,0.37,0.00,0.52,0.27,0.60,0.67,0.27,0.58,0.81,0.41,0.49,0.49,1.0,0.76,0.5])
		return amino_acid_names, hydrophobicities

	elif scale_name == "kd":
						   #Ala,Arg,Asn,Asp,Cys,Gln,Glu,Gly,His,Ile,Leu,Lys,Met,Phe,Pro,Ser,Thr,Trp,Tyr,Val
		hydrophobicities = np.array([0.70,0.0,0.11,0.11,0.78,0.11,0.11,0.46,0.14,1.0,0.92,0.07,0.71,0.81,0.32,0.41,0.42,0.40,0.36,0.97])
		return amino_acid_names, hydrophobicities

	elif scale_name == "eis":
		hydrophobicities = np.array([0.73,0.0,0.40,0.56,0.65,0.39,0.41,0.70,0.49,0.90,0.94,0.24,0.73,0.86,0.61,0.54,0.57,0.77,0.64,1.0])
		return amino_acid_names, hydrophobicities

	elif scale_name == "count":
		"A","I","L","F","W","V"
		hydrophobicities = np.array([1.,0.,0.,0.,0.,0.,0.,0.,0.,1.,1.,0.,0.,1.,0.,0.,0.,1.,0.,1.])
		return amino_acid_names, hydrophobicities



def calculate_hydro_charge(sequence, scale_name='ghavami'):
	aa_names,hydrophobicities = pick_scale(scale_name)

	#aa_index = [aa_names.index(residue) for residue in sequence]

	local_hydrophobicities = [hydrophobicities[index] for index in [aa_names.index(residue) for residue in sequence]]
	indices = [i+1 for i in range(len(sequence))]
	counter_charged = 0
	for aa in sequence:
		if aa == 'D' or aa == 'E' or aa == 'K' or aa == 'R' :
			counter_charged+=1.
	return local_hydrophobicities,index,float(counter_charged),

def calculate_CH_ratio(sequence, scale_name='ghavami'):
	aa_names,hydrophobicities = pick_scale(scale_name)

	#aa_index = [aa_names.index(residue) for residue in sequence]

	local_hydrophobicities = [hydrophobicities[index] for index in [aa_names.index(residue) for residue in sequence]]
	indices = [i+1 for i in range(len(sequence))]
	counter_charged = 0
	for aa in sequence:
		if aa == 'D' or aa == 'E' or aa == 'K' or aa == 'R' :
			counter_charged+=1.
	h_dens = sum(local_hydrophobicities)/len(local_hydrophobicities)
	c_dens = counter_charged/len(local_hydrophobicities)
	ch_ratio = c_dens/h_dens
	return h_dens,c_dens,ch_ratio

def calculate_net_charge(sequence):
	counter_neg=0.
	counter_pos=0.
	for aa in sequence:
		if aa == 'D' or aa == 'E':
			counter_neg+=1.
		if aa == 'K' or aa == 'R':
			counter_pos+=1.
	net_charge = -1*counter_neg+counter_pos
	ncpr = net_charge/len(sequence)
	return net_charge,ncpr


dirname=os.getcwd()

ext = '.txt'



f = open('output_CH_other_NupY.dat','w')

f.write('Name\t\tC/H\tNC\tNCPR\n')

for files in os.listdir(dirname):
	if files.endswith(ext):


#path_file_org = sys.argv[1]
		original_org = getseq(files,1,1)
		collapsed_org = original_org[:612] #EXCLUDES THE FINAL ONE. THIS IS HARDCODED.
		#nupY_extended_org = nupY_original_org[612:]

		full_filepath=dirname+ '/' + files

		prot_name_full=full_filepath.split('NupY-sequences/')[1]
		prot_name=prot_name_full.split('.txt')[0]
		dFG=prot_name.split('_')[1]
		CH =prot_name.split('_')[2]

		hydrophobicities_or,indices_or,charges_or = calculate_hydro_charge(collapsed_org)
		h_dens_or,c_dens_or,ch_ratio_or = calculate_CH_ratio(collapsed_org)
		nc_old,ncpr_old = calculate_net_charge(collapsed_org)

		#h_dens_ext_or, c_dens_ext_or, ch_ratio_ext_or = calculate_CH_ratio(nupY_extended_org)
		f.write('%s\t%s\t%.5f\t%2.0f\t%.5f\n' % (dFG,CH, ch_ratio_or, nc_old, ncpr_old) )

f.close()