#analyze_sequences.py
#(C) 2020 HW de Vries, RUG
# This piece of python code calculates various sequence properties. These include:
# - Percentage of specific types of amino acids within an amino acid sequence
# - Hydrophobicity (as C/H ratio, moving average, varying definitions/scales)
# - Charge patterning, net charge per residue, charge density, Like-charged regions
# - Hydrophobicity patterning score
# - The number and patterning of backbone-stiffness-modifying residues such as Gly and Pro
# 
# Several options can be passed along, such as:
#     -  the division of the protein in different domains
#     -  the comparison between the given sequence and a template (reference) sequence



import sys
import numpy as np
import random
from math import ceil
import copy
import collections
import matplotlib.pyplot as plt


def aa_names():
	return ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

def makestring(seq):
	str_tmp=""
	for x in seq:
		str_tmp+=x
	return str_tmp
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
	def check_charge(res):
		if res == 'D' or res == 'E':
			charge_value = 0-1.

		elif res == 'K' or res == 'R':
			charge_value = 1
		else:
			charge_value = 0

		return charge_value

	local_charges = [ check_charge(residue) for residue in sequence]

	indices = [i+1 for i in range(len(sequence))]
	counter_charged = 0
	for aa in sequence:
		if aa == 'D' or aa == 'E' or aa == 'K' or aa == 'R' :
			counter_charged+=1.
	return local_hydrophobicities, local_charges, index,float(counter_charged)

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



def calculate_AA_histograms(ref_sequence, analysis_sequence, exclude_FG_GLFG = True, normalize_output = True, store_graph = False ):

	'''
	calculate_AA_histograms:
	Inputs: 
	ref_sequence (str): sequence of reference protein domain
	analysis_sequence (str): sequence of analysis protein domain
	exclude_FG_GLFG (boolean, default=True): remove FG/GLFG motifs from the sequence analysis.
	normalize_output (boolean, default=True): return AA vectors that are normalized such that the sum equals 1.
	store_graph (boolean, default=True): store a graph of the histograms, produced by matplotlib. 
										 Note that the output is a normalized histogram by default.!

	Outputs:
	occurrences_reference_protein (np.array)
	occurrences_analysis_protein (np.array)
	... .png (stored in ../..)


	'''
	
	#turn sequences into collections objects:


	reference_protein_counter_obj = collections.Counter(ref_sequence)
	analysis_protein_counter_obj  = collections.Counter(analysis_sequence)
	
	#import aa names
	aa_list = aa_names()


	occurrences_reference_protein=np.array([float(reference_protein_counter_obj[aa]) for aa in aa_list])
	occurrences_analysis_protein=np.array([float(analysis_protein_counter_obj[aa]) for aa in aa_list])
	#for aa in aa_list:
	#	tmp = reference_protein_counter_obj[aa]
#		occurrences_reference_protein.append(tmp)
#		tmp = analysis_protein_counter_obj[aa]
#		occurrences_analysis_protein.append(tmp)
	#use list comprehension to generate np.array with the raw counts of the varying single-letter AA occurrences
	#in the sequences of reference_protein and analysis_protein.


	if normalize_output == True:
		total_AA_ref = float(np.sum(occurrences_reference_protein))
		total_AA_an = float(np.sum(occurrences_analysis_protein))
		occurrences_reference_protein /= total_AA_ref #[occ*/(total_AA_ref) for occ in occurrences_reference_protein]
		occurrences_analysis_protein /= total_AA_an  #[occ*/total_AA_an) for occ in occurrences_analysis_protein]

	#generate matplotlib grap

	fig, ax = plt.subplots(figsize=(10, 6),dpi=300)
	x_range_AAs = np.linspace(1,20,20)
	width=0.35

	rects1 = ax.bar(x_range_AAs - width/2, occurrences_reference_protein*100, width, label=reference_protein_name)
	rects2 = ax.bar(x_range_AAs + width/2, occurrences_analysis_protein*100, width, label=analysis_protein_name)

	plt.tight_layout()
	ax.legend()
	ax.set_xticks(x_range_AAs)
	ax.set_xticklabels(('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'))
	#plt.xticks([], aa_names)
	ax.grid()
	ax.set_ylabel('Frequency (\%)')


	#ax[0].set_xticklabels()
	#ax[1].set_xticklabels()



        #fig1, ax1 = plt.subplots(nrows = 2, ncols = 1, dpi=300)
        #fig1.set_size_inches(10,6)
        #plt.tight_layout()

        #ax1[0].plot(collapsed_indices,hydro_collapsed,label=scale)
        #ax1[0].set_ylabel('Normalized hydrophobicity')
        #ax1[0].set_title(scale)


        #ax1[1].plot(extended_indices,hydro_extended,label=scale)
        #ax1[1].set_ylabel('Normalized hydrophobicity')
        #ax1[1].set_title(scale)
        #ax1[1].set_xlabel('amino acid index')

        plt.savefig('./AA_histogram.png')

	



	return occurrences_reference_protein, occurrences_analysis_protein


def calculate_percentages_AA_groups(occurrence_vector):

	#occurrence vector is formatted in the same order as the aa_list variable.
	aa_list = aa_names()

	hydrophobic_group = ["A","I","L","F","W","V"] #same def. as yamada 
	polar_uncharged_group = ["N","Q","S","T"]
	aromatic_group = ["F","Y","W"]
	charged_group = ["D","E","K","R"]
	charge_pos_group = ["K","R"]
	charge_neg_group = ["D","E"]



	occurrence_groups=np.zeros((6,1)) #same dimensions as the number of groups that's analyzed
	list_groups = ['Hydrophobic', 'Polar q0', 'Aromatic', 'Charged', 'Pos. charged', 'Neg. charged']

	for index1,group in enumerate([hydrophobic_group, polar_uncharged_group, aromatic_group, charged_group, charge_pos_group, charge_neg_group]):
		occ_scal = 0.

		# Loop over the 20 amino acids, and check whether they are present in the considered group.
		# When 
		for index2,aa in enumerate(aa_list):
			if aa in group:
				occ_scal+= occurrence_vector[index2]


		occurrence_groups[index1]=occ_scal


	output_dict = {grp:float(val) for grp,val in zip(list_groups,occurrence_groups)}

	return output_dict


  




#PARSE ARGUMETNS



domain_divisor = int(sys.argv[2]) #assumed to be universal

reference_protein_path = str(sys.argv[1])
reference_protein_name = reference_protein_path.split('.txt')[0]

reference_protein_sequence = getseq(reference_protein_path,1,1)
reference_protein_domain_sequence_1 = reference_protein_sequence[:domain_divisor] #EXCLUDES THE FINAL ONE. THIS IS HARDCODED.
reference_protein_domain_sequence_2 = reference_protein_sequence[domain_divisor:]
reference_protein_string_1 = makestring(reference_protein_domain_sequence_1)

analysis_protein_path = str(sys.argv[3])
analysis_protein_name = analysis_protein_path.split('.txt')[0]
analysis_protein_sequence = getseq(analysis_protein_path,2,1)
analysis_protein_domain_sequence_1 = analysis_protein_sequence[:domain_divisor] #EXCLUDES THE FINAL ONE. THIS IS HARDCODED.
analysis_protein_domain_sequence_2 = analysis_protein_sequence[domain_divisor:]
analysis_protein_string_1 = makestring(analysis_protein_domain_sequence_1)




occurrences_ref, occurrences_an = calculate_AA_histograms(reference_protein_string_1,analysis_protein_string_1)

occurrences_groups_dict_ref = calculate_percentages_AA_groups(occurrences_ref)
occurrences_groups_dict_an = calculate_percentages_AA_groups(occurrences_an)

#output_filename = 'occurrences_groups.csv'
#with open(output_filename,"w") as file:
#		for grp,occ in occurrences_groups_dict_an.items():
#			file.write("%s,\t%s\n" % (grp, str(occ*100) ) )

output_filename2 = 'occurrences_grps_tmp.dat'
with open(output_filename2,"w") as file:
	for pair in occurrences_groups_dict_an.items():
		file.write("%s,\t" % (pair[0])) #,\t %s,\t %s,\t %s,\t %s,\t %s\n
	file.write("\n")
	for pair in occurrences_groups_dict_an.items():
		file.write("%5f,\t" % (100*pair[1]))


output_filename3 = 'occurrence_vector_aa_tmp.dat'

with open(output_filename3,"w") as file:
	for aa in aa_names():
		file.write("%s,\t" % (aa)) #,\t %s,\t %s,\t %s,\t %s,\t %s\n
	file.write("\n")
	for eps in occurrences_an:
		file.write("%2.2f,\t" % (100*eps))

hydrophobicities_an,charges_an,indices_an,charge_total_an = calculate_hydro_charge(analysis_protein_domain_sequence_1)


output_filename4 = 'averages_tmp.dat'
with open(output_filename4,"w") as file:
	headerstr=["eps_sum","eps_PR","|qtot|","NC","NCPR","qptot","qntot"]
	for h in headerstr:
		file.write("%s,\t" % (h)) #,\t %s,\t %s,\t %s,\t %s,\t %s\n
	file.write("\n")
	eps_sum= sum(hydrophobicities_an)
	eps_PR = eps_sum/len(hydrophobicities_an)

	abs_qtot = sum([abs(i) for i in charges_an])
	NC = sum(charges_an)
	NCPR = NC/len(charges_an)
	qptot= sum([i for i in charges_an if i > 0])
	qntot = sum([-i for i in charges_an if i < 0])

	print_vals=[eps_sum, eps_PR, abs_qtot, NC, NCPR, qptot, qntot]
	for vals in print_vals:
		file.write("%f,\t" % (vals))



#DEBUGGING
#counter=0
#for eps,q in zip(hydrophobicities_or,charges_or): #
#	print(reference_protein_domain_sequence_1[counter],counter+1,eps,q)
#	counter+=1



#calculate C/H ratio
#hydrophobicities_or,indices_or,charges_or = calculate_hydro_charge(REFERENCE_HERE)
#h_dens,c_dens,ch_ratio = calculate_CH_ratio(REFERENCE_HERE)
#nc_old,ncpr_old = calculate_net_charge(REFERENCE_HERE)

#h_dens_ext_or, c_dens_ext_or, ch_ratio_ext_or = calculate_CH_ratio(REFERENCE_HERE)

#print('CH ratio of collapsed domain is: ',ch_ratio_or)
#print('CH ratio of extended domain is: ',ch_ratio_ext_or)

