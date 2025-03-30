import MDAnalysis as mda
import os
import sys
import numpy as np
import random as rd

prot1name=str(sys.argv[1])
prot2name=str(sys.argv[2])
#atomgr_overlap = u2.select_atoms('within 5 of 

#parse_inputs

box_size_L       = int(sys.argv[3])
master_struct  = mda.Universe(prot1name)
ref_struct = mda.Universe(prot2name)



def randpm():
	return 1 if rd.random() < 0.5 else -1

def generate_trial_conf(original_struct, box_size=30, offset_boundary = 4):

	additional_structure=mda.Universe(original_struct.filename)
	#offset boundary set to 4.5 conform radius of gyration.
	gen_x = 10*(rd.random()* (box_size - offset_boundary/2.0) + offset_boundary/2.0)
	gen_y = 10*(rd.random()* (box_size - offset_boundary/2.0) + offset_boundary/2.0)
	gen_z = 10*(rd.random()* (box_size - offset_boundary/2.0) + offset_boundary/2.0)

	new_cent = np.array([gen_x,gen_y,gen_z])
	com = additional_structure.atoms.center_of_geometry()
	diff = com-new_cent

	additional_structure.atoms.positions = additional_structure.atoms.positions - diff
	new_com = additional_structure.atoms.center_of_geometry()

	return additional_structure.atoms

def overlap_check(existing_struct, new_struct,exclusionRadius=6):

	points_within_radius = 0

	#unpack positions of new structure
	for pos in new_struct.atoms.positions:
		atoms_overlap = existing_struct.select_atoms('point '+str(pos[0])+' '+str(pos[1])+' '+str(pos[2])+' '+str(exclusionRadius))
		if len(atoms_overlap) !=0: #the site is too close to existing beads, so the point is occupied -> return True
			points_within_radius+=1 

	if points_within_radius ==0:
		return False
	else: 								#the site is not too close to existing beads, so the point is not occupied -> return False
		return True

def check_OOB(new_struct,box_size=30):
	fcheckpos = new_struct.atoms.positions[:]
	return (fcheckpos > 0).all() and (fcheckpos < 10*box_size).all()
		
num_trials = 0
done_adding=False
while done_adding==False:

	new_conf = generate_trial_conf(ref_struct,box_size_L)
	overlap_bool = overlap_check(master_struct,new_conf)
	within_bounds = check_OOB(new_conf,box_size_L)
	if overlap_bool == False and within_bounds == True : 
		master_struct = mda.Merge(master_struct.atoms,new_conf.atoms)
		done_adding=True
	else:
		num_trials+=1
		print('overlap was:',overlap_bool,'within_bounds was:',within_bounds)




master_struct.atoms.write('yFGKap-A_1_complex.pdb')
os.system('gmx editconf -f yFGKap-A_1_complex.pdb -box 30 30 30 -bt tric -c -o yFGKap-A_1_complex.pdb')

