#!/bin/bash


prot_name=$1
len_prot_1=$2
len_prot_2=$3
box_size=$4

input_struct=${prot_name}-A_1_complex.pdb




#0. gen ndx groups
python gen_groups.py $input_struct $len_prot_1 $len_prot_2

#1. Do polystat to determine radii of gyration. We will define the 'subvolume' as the interacting volume. That is:
#   the volume with a diameter equal to the average of the to largest diameters of the proteins.
#   I.O.W.: D1 = 2xRg(1), D2 = 2xRg(2), R_subv = (D1 + D2)/2

cp replica_1/MD.tpr .
cp replica_1/yFGKap-A_1_complex.pdb .

echo "PROT1" | gmx polystat -f traj_combined.xtc -s MD.tpr -n kd.ndx -o prot_1_polystat.xvg -xvg none
echo "PROT2" | gmx polystat -f traj_combined.xtc -s MD.tpr -n kd.ndx -o prot_2_polystat.xvg -xvg none


#2. Calculate, using gmx pairdist, the minimum distances between all the atoms, and the distances between the COM of the proteins.

#just the -minimum- distance between any of the residues
gmx pairdist -f traj_combined.xtc -n kd.ndx -s MD.tpr -ref PROT1 -sel PROT2 -pbc -type min -o pairdist_min.xvg -xvg none
#the COM distance betwen both proteins
gmx pairdist -f traj_combined.xtc -s MD.tpr -n kd.ndx -ref PROT1 -sel PROT2 -selrpos whole_mol_com -seltype whole_mol_com -o pairdist_com.xvg -xvg none

#the max internal distance, used for averaging.
gmx pairdist -f traj_combined.xtc -n kd.ndx -s MD.tpr -ref PROT2 -sel PROT2 -pbc -type max -o pairdist_prot_2_max.xvg -xvg none

rm kd.dat
python calculate_kd_alt_idp.py $box_size $len_prot_1 $len_prot_2 $input_struct
