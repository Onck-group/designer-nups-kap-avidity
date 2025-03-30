#!/bin/bash

for val in {1..10} ; do

	cd replica_${val}
	cp ../js_PG1 .
	cp ../ext_PG1 .
	sbatch js_PG1 yFGKap-A_1_complex.pdb
	cd ..
done 

for val in {11..20} ; do
	cd replica_${val}
	cp ../js_PG2 .
	cp ../ext_PG2 .
	sbatch js_PG2 yFGKap-A_1_complex.pdb
	cd ..
done
