#!/bin/bash

for dFG in 7 13 26 52 104 ; do
	
	for ch_rat_dir in reference_structures/${dFG}/* ; do
		
		ch_val=$( echo $ch_rat_dir | cut -f 3 -d "_" )
		
		mkdir -p NupY_${dFG}_${ch_val}
		cp -r input_files/* NupY_${dFG}_${ch_val}
		cp $ch_rat_dir/NupY*.pdb NupY_${dFG}_${ch_val}/basefiles/
		cd NupY_${dFG}_${ch_val}/basefiles/
		python rebuild_Nup_G1.py NupY_${ch_val}*.pdb NupY_${dFG}_${ch_val} True True
		cat NupY_${dFG}_${ch_val}_rebuilt.pdb > yFGKap-B.pdb
		cat NupY_${dFG}_${ch_val}.itp > yFGKap-B.itp
		sed -i "s/NupY_${dFG}_${ch_val}/yFGKap-B/g" yFGKap-B.itp	
		python combine_mols.py yFGKap-A.pdb yFGKap-B.pdb 45
	        python gen_groups_NupY.py yFGKap-A_1_complex.pdb 861 803	
		rm *\#
		cd ..
		sed -i "s/XX.XX/${dFG}_${ch_val}/g" ext_PG1
		sed -i "s/XX.XX/${dFG}_${ch_val}/g" ext_PG2
		sed -i "s/XX.XX/${dFG}_${ch_val}/g" js_PG1
		sed -i "s/XX.XX/${dFG}_${ch_val}/g" js_PG2
		cd ..
	done

	

done
