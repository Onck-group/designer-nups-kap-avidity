import MDAnalysis as mda
import sys

input_file = sys.argv[1]
len_prot_1 = int(sys.argv[2])
len_prot_2 = int(sys.argv[3])
u = mda.Universe(input_file)

if len(u.atoms) != len_prot_1+len_prot_2:
	print('WARNING: lengths not equal!')
	print('will exit the gen_groups.py script now. Please check your numbering!')
	exit()

with mda.selections.gromacs.SelectionWriter('kd.ndx',mode='w') as ndx:
		ndx.write(u.select_atoms('bynum 1:%s' % (str(len_prot_1))), name = 'PROT1')
		ndx.write(u.select_atoms('bynum %s:%s' %  (str(len_prot_1+1), str(len_prot_1+len_prot_2) )), name = 'PROT2')
