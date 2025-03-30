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

with mda.selections.gromacs.SelectionWriter('sys.ndx',mode='w') as ndx:
		ndx.write(u.select_atoms('all'), name = 'SYSTEM')
		ndx.write(u.select_atoms('bynum 1:%s' % (str(len_prot_1))), name = 'PROT1')
		ndx.write(u.select_atoms('bynum %s:%s' %  (str(len_prot_1+1), str(len_prot_1+len_prot_2) )), name = 'PROT2')
		ndx.write(u.select_atoms('name B*'), name = 'BS')
		ndx.write(u.select_atoms('name F F1'), name = 'F')
		ndx.write(u.select_atoms('name F1'), name = 'F1')
		ndx.write(u.select_atoms('name G1'), name = 'G1')
		ndx.write(u.select_atoms('bynum 1:%s' % (str(len_prot_1))), name = 'FREEZ')
                
		#cation-pi interactions
		ndx.write(u.select_atoms('name R RR K KK'),name='RK')
                ndx.write(u.select_atoms('name F FF Y YY W WW'), name = 'FYW')

                ndx.write(u.select_atoms('name BR BK'),name='BRK')
                ndx.write(u.select_atoms('name BF BY BW'),name='BFYW')

                ndx.write(u.select_atoms('not (name BF F F1 G1 FF BY Y YY BW W WW BR R RR BK K KK)'), name = 'AA')





