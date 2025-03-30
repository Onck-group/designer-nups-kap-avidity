import numpy as np 
import MDAnalysis as mda
import MDAnalysis.analysis.distances
import os
import sys






polystat_filenames = ["prot_1_polystat.xvg", "prot_2_polystat.xvg"]
mindist_filename = "pairdist_min.xvg"
com_filename     = "pairdist_com.xvg"

box_size = int(sys.argv[1])
prot_1_size = int(sys.argv[2])
prot_2_size = int(sys.argv[3])
input_struct = sys.argv[4]





def read_file(file_name,column_no):

		if column_no == 'None':
			arr_output = [line.split() for line in open(file_name)]
		else:
			arr = [line.split() for line in open(file_name)]
			arr_output = np.array( [arr[timestep][column_no] for timestep in range(len(arr))] ).astype(np.float)

		return arr_output


def calc_rg(filename_list):
	'''
	calc_rg(filename_list): This subroutine calculates the Rg and associated standard deviation for protein_1 and protein_2.
	This assumes that you prepared two polystat.xvg files without headers.
	POlystat files are indexed as:

	[0]       [1]    [2]    [3,4,5]
	timestep  D_ee    Rg    PC of Rg tensor
	'''




	rg_values = np.zeros(len(filename_list))
	rg_std    = np.zeros(len(filename_list))




	for index,file_name in enumerate(filename_list):

		#layout of polystat.xvg is typically: time / D_ee / Rg, so we use column three, which is no. 2 in zero-based indexing.
		rg_array = read_file(file_name,2)
   
		mean_rg = np.mean(rg_array)
		std_value = np.std(rg_array)
		rg_values[index] = mean_rg
		rg_std[index] = std_value


	return rg_values, rg_std



def calc_bound_fraction(cut_off_value, input_file, method_type = 'mindist'):
	'''
	Calculates: 
	- The number of unbound states N0
	- The number of bound states N1
	- pb, which is N1/N_tot
	'''

	if method_type == 'mindist':
		minimum_distance_array = read_file(input_file,column_no = 1)

		tot_frames = len(minimum_distance_array)

		N_unbound = 0.
		N_bound = 0.
		for dist in minimum_distance_array:
			if dist <= cut_off_value:
				N_bound +=1
			else:
				N_unbound +=1

		#check for consistency
		N_tot = N_unbound+N_bound
		if N_tot != tot_frames:
			print('WARNING! Tot. frames is not equal to sum of fractions')
			exit()

		bound_fraction = N_bound/N_tot





	return N_unbound,N_bound,bound_fraction



def calc_subvolume_radius(protein_filename, prot_length_1, prot_length_2):

	u = mda.Universe(protein_filename)
	prot_1 = u.select_atoms('bynum 1:%s' % str(prot_length_1))
	prot_2 = u.select_atoms('bynum %s:%s' % ( str(prot_length_1+1), str(prot_length_1+prot_length_2) ))
	if len(prot_1)+len(prot_2) != len(u.atoms) :
		print('WARNING: issue with protein lengths, respectively')
	d_prot_1 = MDAnalysis.analysis.distances.self_distance_array(prot_1.positions)
	d_prot_2 = MDAnalysis.analysis.distances.self_distance_array(prot_2.positions)

	D_1 = np.max(d_prot_1)
	D_2 = np.max(d_prot_1)

	#Definition as given by Lopez et al. JPhysChemB 2020: R_sub = avg(D1,D2) with D1,D2 the largest diameters of the proteins
	subvolume_radius = (D_1+D_2)/(2*10) 
	# *10 is included to convert Angstrom output of MDA to nm --> GMX compatibility


	return subvolume_radius


def calc_subvolume_fraction(cut_off_value, input_file):
	'''
	'''


	distance_array = read_file(input_file,column_no = 1)

	total_frames = len(distance_array)

	N_outside = 0.
	N_inside = 0.
	for dist in distance_array:
		if dist <= cut_off_value:
			N_inside +=1
		else:
			N_outside +=1

	N_total = N_inside+N_outside
	if N_total != total_frames:
		print('WARNING! Tot. frames is not equal to sum of fractions')
		exit()


	subvolume_fraction = N_inside / N_total

	return subvolume_fraction


def calc_kd(box_volume, dimerization_volume, subvolume, N_unbound, N_bound, subvolume_fraction, bound_fraction, method_type = 'Lopez'):
	'''
	
	'''
	
	N_A = 6.02214076*10**(23) # mol^-1
	conversion_to_molar = 10**(24) #from nm^3 to L


	if method_type == 'Kim':
		k_d = conversion_to_molar * ( (1 - bound_fraction )**2 / (N_A * box_volume * bound_fraction) )
		print(k_d)
	if method_type == 'DeJong':
		k_d = conversion_to_molar * ( N_unbound / (N_bound * N_A * (box_volume - dimerization_volume) ) )
		print(k_d)
	if method_type == 'Lopez':
		k_d = conversion_to_molar * ( 1/( N_A * box_volume * bound_fraction ) ) * ( (1-subvolume_fraction) / (1-subvolume/box_volume) )
		print(k_d)
	return k_d





rg_values,rg_std = calc_rg(filename_list = polystat_filenames)
#rg_upper_lim = [rg+std for rg,std in zip(rg_values,rg_std)]
#rg_lower_lim = [rg-std for rg_std in zip(rg_values,rg_std)]

R_subvolume = calc_subvolume_radius(input_struct,prot_1_size,prot_2_size)

#OLD: (2*rg_values[0]+2*rg_values[1])/2
#R_subvolume_upper = (2*rg_upper_lim[0]+2*rg_upper_lim[1])/2
#R_subvolume_lower = (2*rg_lower_lim[0]+2*rg_lower_lim[1])/2

V_subvolume = 4/3*np.pi*R_subvolume**3 #in nm^3
#V_subvolume_upper = 4/3*np.pi*R_subvolume_upper**3
#V_subvolume_lower = 4/3*np.pi*R_subvolume_lower**3

V_box = box_size**3 #in nm^3
V_dimerization = 0.

subvolume_fraction = calc_subvolume_fraction(R_subvolume, com_filename)
#subvolume_fraction_upper = calc_subvolume_fraction(R_subvolume_upper, com_filename)
#subvolume_fraction_lower = calc_subvolume_fraction(R_subvolume_lower, com_filename)

N_unbound, N_bound, bound_fraction = calc_bound_fraction(cut_off_value = 0.8, input_file = mindist_filename, method_type= 'mindist')

kd = calc_kd(box_volume = V_box, dimerization_volume = V_dimerization, subvolume = V_subvolume, N_unbound = N_unbound, N_bound = N_bound, subvolume_fraction = subvolume_fraction, bound_fraction = bound_fraction, method_type = 'Lopez')
#kd_upper = calc_kd(box_volume = V_box, dimerization_volume = V_dimerization, subvolume = V_subvolume_upper, N_unbound = N_unbound, N_bound = N_bound, subvolume_fraction = subvolume_fraction_upper, bound_fraction = bound_fraction, method_type = 'Lopez')
#kd_lower = calc_kd(box_volume = V_box, dimerization_volume = V_dimerization, subvolume = V_subvolume_lower, N_unbound = N_unbound, N_bound = N_bound, subvolume_fraction = subvolume_fraction_lower, bound_fraction = bound_fraction, method_type = 'Lopez')



curr_eps = float(os.getcwd().split('syst_')[1])

f = open('kd.dat','w+')


f.write('%.3f\t%f\t%f\t%f\t%f\n' % (curr_eps, kd*10**6, bound_fraction, subvolume_fraction,V_subvolume/V_box))
f.close()





