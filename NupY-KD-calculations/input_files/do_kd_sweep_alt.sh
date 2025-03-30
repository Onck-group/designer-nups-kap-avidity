#!/bin/bash


prot_name=$1
len_prot_1=$2
len_prot_2=$3
box_size=$4

input_struct=${prot_name}-A_1_complex.pdb

rm kd_output_pb_pv_alt.dat

for simdir in NupY_* ; do
	cd $simdir
	cp replica_1/MD.tpr .
	cp ../analysis_files/* .
	
	sh do_kd_calculation_alt.sh $prot_name $len_prot_1 $len_prot_2 $box_size 
	cd ..
	
	cat $simdir/kd.dat >> kd_output_pb_pv_alt.dat	
done

