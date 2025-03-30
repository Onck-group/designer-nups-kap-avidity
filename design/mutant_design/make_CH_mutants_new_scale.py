import sys
import numpy as np
import random
from math import ceil
import copy



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
		
		hydrophobicities = np.array([0.7,0.05,0.33,0.05,0.68,0.64,0.05,0.41,0.53,0.98,1,0.05,0.78,1.00,0.65,0.45,0.51,0.96,0.82,0.94])
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


def replace_FG_GLFG_motif(input_sequence,toggle_aromatic=True):
	'''
	Takes as an input the sequence, finds the FG/GLFG motifs and replaces them with a combination of
	amino acids that have the same overall C/H ratio. The function works in the
	following way:
	- it determines te FG-type
	- then chooses a pair of aas that is a good replacement, based on:
		- hydrophobicity MATCH
		- the number of aromatic residues
	- two 'pairs' are chosen to that end

	'''

	aromatic_list= [['W','A'],['W','C']]
	alternative_list=[['C','V'],['I','C'],['A','V']]
	output_sequence = input_sequence[:]
	#BUG: INSERTION DOES NOT WORK CORRECTLY
	for index,residue in enumerate(input_sequence):

		glfg_selector = (residue == 'F' and input_sequence[index+1] == 'G' and input_sequence[index-1] == 'L' and input_sequence[index-2] == 'G')
		fg_selector = (residue == 'F' and input_sequence[index+1] == 'G' and input_sequence[index-1] != 'L')
		if glfg_selector:
			if toggle_aromatic:
				replacement_1 = random.choice(alternative_list)
				replacement_2 = random.choice(aromatic_list)
				output_sequence[index-2:index]	= replacement_1
				output_sequence[index:index+2]  	= replacement_2

			else:
				replacement_1 = random.choice(alternative_list)
				replacement_2 = random.choice(alternative_list)
				output_sequence[index-2:index]	= replacement_1
				output_sequence[index:index+2]  	= replacement_2

		if fg_selector:
			if toggle_aromatic:
				replacement = random.choice(aromatic_list)
				output_sequence[index:index+2] = replacement
			else:
				replacement = random.choice(alternative_list)
				output_sequence[index:index++2] = replacement
	return output_sequence

def make_CH_spacing_mutant(input_sequence, target_ch_ratio, target_nc, tolerance=0.0005, scale_name = 'ghavami' ):
	'''
	takes as an input sequence and makes random mutations to amino acids such that a designated CH ratio is obtained
	'''

	#set several conditions: the default set of AAs is a set of disorder promoting/no preference AAs with no charged residues.
	#we prefer to modulate the hydrophobicity and not so much the charge density (which is already low).

	aa_names,hydrophobicities = pick_scale(scale_name)
	set_of_AAs = ['A','G','Q','S','P','M','T','H','E','K','R','D']
	set_of_hydrophobicities = [hydrophobicities[index] for index in [aa_names.index(residue) for residue in set_of_AAs]]
	#alternative_set_of_AAs:
	output_sequence = input_sequence[:]


	current_h_dens, current_c_dens, current_ch_ratio = calculate_CH_ratio(output_sequence) #output_sequence is iteratively updated.
	epsilon = (current_ch_ratio - target_ch_ratio)
	offset_sign = np.sign(epsilon)
	print('initializing CH ratio mutation:')
	print('current CH ratio:', current_ch_ratio)
	print('target CH ratio:',target_ch_ratio)
	print('initial epsilon:', epsilon)

	indices_FG_or_GLFG=[]
	for index,residue in enumerate(input_sequence):

		glfg_selector = (residue == 'F' and input_sequence[index+1] == 'G' and input_sequence[index-1] == 'L' and input_sequence[index-2] == 'G')
		fg_selector = (residue == 'F' and input_sequence[index+1] == 'G' and input_sequence[index-1] != 'L')   
		
		if glfg_selector:
			placeholder =range(index-2,index+2)
			for i in placeholder:
				indices_FG_or_GLFG.append(i)

		if fg_selector:

			placeholder =range(index,index+2)
			for i in placeholder:
				indices_FG_or_GLFG.append(i)
	spacer_indices = [i for i in range(len(input_sequence)) if i not in indices_FG_or_GLFG]

	while abs(epsilon) > tolerance :
		random_index = random.choice(spacer_indices)
		random_residue = output_sequence[random_index]
		print('random residue is:', random_residue)
		tmp_index = aa_names.index(random_residue)

		random_residue_hydrophobicity = hydrophobicities[tmp_index] 

		if offset_sign > 0 :
			available_residues = [residue for index,residue in enumerate(set_of_AAs) if set_of_hydrophobicities[index] > random_residue_hydrophobicity]
			if len(available_residues) > 0:
				output_sequence[random_index] = random.choice(available_residues)
				print('turned residue ',random_residue, 'into:', output_sequence[random_index])

		if offset_sign < 0 :
			available_residues = [residue for index,residue in enumerate(set_of_AAs) if set_of_hydrophobicities[index] < random_residue_hydrophobicity]
			print(available_residues)
			if len(available_residues) > 0:
				new_residue = random.choice(available_residues)
				nc_before,ncpr_before = calculate_net_charge(output_sequence)
				offset_nc_before = (nc_before - target_nc)/target_nc
				output_sequence[random_index] = new_residue
				nc_after,ncpr_after = calculate_net_charge(output_sequence)
				offset_nc_after = (nc_after - target_nc)/target_nc
				if abs(offset_nc_after) <= abs(offset_nc_before):
					pass
				else:
					available_residues = [item for item in available_residues if item is not new_residue]
					if len(available_residues) > 0:
						new_residue = random.choice(available_residues)
						output_sequence[random_index] = new_residue

				print('turned residue ',random_residue, 'into:', output_sequence[random_index])

		current_h_dens, current_c_dens, current_ch_ratio = calculate_CH_ratio(output_sequence)
		epsilon = (current_ch_ratio - target_ch_ratio)
		offset_sign = np.sign(epsilon)
		print('Epsilon currently is:', epsilon)

	return output_sequence



def make_FG_spacing_mutant(input_sequence, new_spacing):
	domain_length = len(input_sequence)
	dummy_domain = input_sequence[:]
	replaced_iso_hydrophobic_domain = replace_FG_GLFG_motif(input_sequence)
	output_domain = replaced_iso_hydrophobic_domain[:]
	indices_original_FG=[]
	indices_original_GLFG=[]
	print('working on loop')
	for index,residue in enumerate(input_sequence):

		glfg_selector = (residue == 'F' and input_sequence[index+1] == 'G' and input_sequence[index-1] == 'L' and input_sequence[index-2] == 'G')
		fg_selector = (residue == 'F' and input_sequence[index+1] == 'G' and input_sequence[index-1] != 'L')
		if glfg_selector:
			dummy_domain[index-2:index+2]=['DUM','DUM','DUM','DUM']

			indices_original_GLFG.append(range(index-2,index+2))

		if fg_selector:
			dummy_domain[index:index+2]=['DUM','DUM']
			indices_original_FG.append(range(index,index+2))


	#determine the motifs for the spacers.
	FG_indices_new_spacing, GLFG_indices_new_spacing = determine_indices_motifs(domain_length, new_spacing)
	replacable_FG_indices = indices_original_FG[:]
	replacable_GLFG_indices = indices_original_GLFG[:]

	conflict_indices=[]
	re_insert_motifs_FG=[]
	placeholder=[]
	for number,index in enumerate(FG_indices_new_spacing):

		for pair in indices_original_FG:
			if index in pair:
				conflict_indices.append(index)
				if number % 2 ==0:
					print('conflict in residue ',index, 'residue: ', replaced_iso_hydrophobic_domain[index], 'new residue: F',)
				else:
					print('conflict in residue ',index, 'residue: ', replaced_iso_hydrophobic_domain[index], 'new residue: G',)
		if number % 2 == 0:
			output_domain[index] = 'F'
			placeholder.append(replaced_iso_hydrophobic_domain[index])
		else:
			output_domain[index] = 'G'
			placeholder.append(replaced_iso_hydrophobic_domain[index])
		if len(placeholder)==2:
			re_insert_motifs_FG.append(placeholder)
			placeholder=[]


	re_insert_motifs_GLFG=[]
	placeholder=[]
	for number,index in enumerate(GLFG_indices_new_spacing):

		for pair in indices_original_GLFG:
			if index in pair:
				conflict_indices.append(index)
				if number % 4 == 0 :
					print('conflict in residue ',index, 'residue: ', replaced_iso_hydrophobic_domain[index], 'new residue: G',)
				if number % 4 == 1 :
					print('conflict in residue ',index, 'residue: ', replaced_iso_hydrophobic_domain[index], 'new residue: L',)
				if number % 4 == 2 :
					print('conflict in residue ',index, 'residue: ', replaced_iso_hydrophobic_domain[index], 'new residue: F',)
				if number % 3 == 0 :
					print('conflict in residue ',index, 'residue: ', replaced_iso_hydrophobic_domain[index], 'new residue: G',)
		if number % 4 == 0 :
			output_domain[index] = 'G'
			placeholder.append(replaced_iso_hydrophobic_domain[index])
		if number % 4 == 1 :
			output_domain[index] = 'L'
			placeholder.append(replaced_iso_hydrophobic_domain[index])
		if number % 4 == 2 :
			output_domain[index] = 'F'
			placeholder.append(replaced_iso_hydrophobic_domain[index])
		if number % 4 == 3 :
			output_domain[index] = 'G'
			placeholder.append(replaced_iso_hydrophobic_domain[index])
		if len(placeholder) == 4:
			re_insert_motifs_GLFG.append(placeholder)
			placeholder=[]

	replacable_FG_indices = [FG_pair for FG_pair in replacable_FG_indices if not any(item in FG_pair for item in conflict_indices)]
	random.shuffle(replacable_FG_indices)

	replacable_GLFG_indices = [GLFG_pair for GLFG_pair in replacable_GLFG_indices if not any(item in GLFG_pair for item in conflict_indices)]
	random.shuffle(replacable_GLFG_indices)
	print(re_insert_motifs_FG)

	print(re_insert_motifs_GLFG)
	for index,item in enumerate(re_insert_motifs_FG):
		placement_indices = replacable_FG_indices[index]
		print('placement_indices are:',placement_indices)
		for number,residue in enumerate(item):
			print('replacing residue',residue,'at location',placement_indices[number])
			output_domain[placement_indices[number]] = residue
	for index,item in enumerate(re_insert_motifs_GLFG):
		placement_indices = replacable_GLFG_indices[index]
		print('placement_indices are:',placement_indices)
		for number,residue in enumerate(item):
			print('replacing residue',residue,'at location',placement_indices[number])
			output_domain[placement_indices[number]] = residue

	return output_domain


	#

def determine_indices_motifs(domain_length, spacing):

	FG_repeat_length=2
	GLFG_repeat_length=4
	spacer_length=spacing
	repeat_length=4*FG_repeat_length+3*GLFG_repeat_length+7*spacer_length
	#_FG_GLFG_FG_GLFG_FG_GLFG_FG
	#1  2    3  4    5  6    7
	fg_counter=0
	glfg_counter=0
	repeat_counter=0

	#below, we account for zero-based indexing: slot after spacer_length is actually spacer_length in 0-based indexing
	FG_1_start = spacer_length
	FG_1_indices = range(FG_1_start,FG_1_start+2)

	GLFG_1_start = FG_1_start+FG_repeat_length + spacer_length
	GLFG_1_indices = range(GLFG_1_start,GLFG_1_start+4)

	FG_2_start = GLFG_1_start + GLFG_repeat_length + spacer_length
	FG_2_indices = range(FG_2_start,FG_2_start+2)

	GLFG_2_start = FG_2_start+FG_repeat_length + spacer_length
	GLFG_2_indices = range(GLFG_2_start,GLFG_2_start+4)

	FG_3_start = GLFG_2_start + GLFG_repeat_length + spacer_length
	FG_3_indices = range(FG_3_start,FG_3_start+2)

	GLFG_3_start = FG_3_start+FG_repeat_length + spacer_length
	GLFG_3_indices = range(GLFG_3_start,GLFG_3_start+4)

	FG_4_start = GLFG_3_start + GLFG_repeat_length + spacer_length
	FG_4_indices = range(FG_4_start,FG_4_start+2)

	FG_indices_single_repeat = [FG_1_indices, FG_2_indices, FG_3_indices, FG_4_indices]
	GLFG_indices_single_repeat = [GLFG_1_indices, GLFG_2_indices, GLFG_3_indices]
	num_repeats = int(ceil(domain_length/repeat_length))

	FG_indices_full_sequence = []
	GLFG_indices_full_sequence = []

	for i in range(num_repeats+1):
		for motifs in FG_indices_single_repeat:
			print(motifs)
			for index in motifs:
				index += i*repeat_length
				if index <= domain_length-1:
					FG_indices_full_sequence.append(index)

		for motifs in GLFG_indices_single_repeat:

			for index in motifs:
				index+= i*repeat_length
				if index <= domain_length-1:
					GLFG_indices_full_sequence.append(index)


	print('new FG indices:', FG_indices_full_sequence)
	print('new GLFG_indices:', GLFG_indices_full_sequence)
	return FG_indices_full_sequence,GLFG_indices_full_sequence




	hydro_dens_collapsed     = sum(hydro_collapsed)/(len(hydro_collapsed))
	charge_dens_collapsed    = calculate_charge_count(nupY_collapsed)/len(hydro_collapsed)
	ch_ratio_collapsed = charge_dens_collapsed/hydro_dens_collapsed



#MAIN CODE ##############################################
path_file = sys.argv[1]

nupY_original = getseq(path_file,1,1)
nupY_collapsed = nupY_original[:612] #EXCLUDES THE FINAL ONE. THIS IS HARDCODED.
nupY_extended = nupY_original[612:]



#num_variations = 25

#spacing = 52

hydrophobicities_or,indices_or,charges_or = calculate_hydro_charge(nupY_collapsed)
h_dens_or,c_dens_or,ch_ratio_or = calculate_CH_ratio(nupY_collapsed)
nc_old,ncpr_old = calculate_net_charge(nupY_collapsed)

#targets=[0.0239] # 0.0956, 0.1912, 0.2868, 0.3824, 0.478, 0.5736]
targets = [0.0239, 0.1912, 0.3824, 0.5736]#0.0956]
variations = 10000

for target in targets:
	counter = 0
	while counter <= variations:

		nupY_collapsed_new_CH = make_CH_spacing_mutant(nupY_collapsed, target,nc_old)
		h_dens_new,c_dens_new,ch_ratio_new = calculate_CH_ratio(nupY_collapsed_new_CH)
		nc_new,ncpr_new = calculate_net_charge(nupY_collapsed_new_CH)

		diff_nc = (nc_new-nc_old)*100/nc_old

		if diff_nc <= 2.5:
			counter+=1
			NupY_mutant = ["G"] + ["P"] + nupY_collapsed_new_CH+nupY_extended + ["C"]

			prot_name='outputs/1bpa_seqs/104_spacing/NupY_mutant_ch'+str(target)+'var'+str(counter)+'.txt'

			with open(prot_name,"w") as file:
				file.write("%%NupY_ch%f_%d\n" % (target,counter)) 
				for residue in NupY_mutant:
					file.write("%s" % residue)


			prot_name='outputs/fasta/104_spacing/NupY_mutant_ch'+str(target)+'var'+str(counter)+'.fasta'

			with open(prot_name,"w") as file:
				file.write("> NupY_ch%f_%d\n" % (target,counter))
				counter2 = 0.
				for residue in NupY_mutant:
					file.write("%s" % residue)
					counter2+=1
					if counter2 % 60 == 0:
						file.write("\n")












#counter = 0
#for residue1,residue2 in zip(nupY_collapsed,nupY_collapsed_new_CH):
#	print(residue1,residue2,counter)
#	counter+=1


# for spacing in spacings:
# 	for variation in range(1,numVariations+1):
# 		new_NupY_c = make_FG_spacing_mutant(nupY_collapsed, spacing)
# 		hydrophobicities_new,indices_new,charges_new = calculate_hydro_charge(new_NupY_c)

# 		h_dens_or = sum(hydrophobicities_or)/len(hydrophobicities_or)
# 		c_dens_or = charges_or/len(hydrophobicities_or)
# 		ch_ratio_or = c_dens_or/h_dens_or

# 		h_dens_new = sum(hydrophobicities_new)/len(hydrophobicities_new)
# 		c_dens_new = charges_new/len(hydrophobicities_new)

# 		ch_ratio_new = c_dens_new/h_dens_new

# 		differences_or_new = (ch_ratio_new - ch_ratio_or) / ch_ratio_or

# 		print('Old hydrophobicity:',h_dens_or, 'Old charge count:', c_dens_or)
# 		print('New hydrophobicity:',h_dens_new, 'New charge count:', c_dens_new)

# 		print('C/H ratio old:', ch_ratio_or, 'C/H ratio new:', ch_ratio_new, 'difference, percentagewise:', differences_or_new, '%')


# 		NupY_mutant = new_NupY_c+nupY_extended

# 		prot_name='outputs/NupY_FG_spacing_mutant_'+str(spacing)+'var_'+str(variation) + '.txt'

# 		with open(prot_name,"w") as file:
# 			file.write("%%NupY_%d_%d\n" % (spacing,variation))
# 			file.write("GP") 
# 			for residue in NupY_mutant:
# 				file.write("%s" % residue)
# 			file.write("C")



#NupYC_test = replace_FG_GLFG_motif(nupY_collapsed)
#counter = 0
#for residue1,residue2 in zip(nupY_collapsed,new_NupY_c):#
#	print(residue1,residue2,counter)
#	counter+=1
