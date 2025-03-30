#!/bin/bash

pattern=$1


#cleanup

rm averages_tmp.dat
rm averages_${pattern}_out.dat
rm occurrences_grps_tmp.dat
rm occurrences_grps_${pattern}_out.dat
rm occurrence_vector_aa_tmp.dat
rm occurrence_vector_aa_${pattern}_out.dat

for f in NupY_$pattern* ; do
	
	python analyze_sequences.py template 612 $f
	tail -n 1 -q averages_tmp.dat >> averages_${pattern}_out.dat
	echo "" >> averages_${pattern}_out.dat
	
	tail -n 1 -q occurrences_grps_tmp.dat >> occurrences_grps_${pattern}_out.dat
	echo "" >> occurrences_grps_${pattern}_out.dat

	tail -n 1 -q occurrence_vector_aa_tmp.dat >> occurrence_vector_aa_${pattern}_out.dat
	echo "" >> occurrence_vector_aa_${pattern}_out.dat

done
