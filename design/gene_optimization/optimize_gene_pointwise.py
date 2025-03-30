import os
import numpy as np
import sys
import random as rd
from copy import deepcopy
import math
import difflib

path_file = sys.argv[1]
print(path_file)





codon_table={ 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    } 

#in the below table, everything below 30% is left out. if not, then the highest % is taken
AA_table = {
	    'I':{'ATC','ATT'}, 
	    'M':{'ATG'},
	    'T':{'ACC'} ,
        'N':{'AAC','AAT'},
        'K':{'AAA'},
        'S':{'AGC'},
        'R':{'CGC','CGT'},
        'L':{'CTG'},
        'P':{'CCG'},
        'H':{'CAC','CAT'},
        'Q':{'CAA','CAG'},
        'V':{'GTG'},
        'A':{'GCG'},
        'D':{'GAC','GAT'},        
        'E':{'GAA','GAG'},
        'G':{'GGC','GGT'},
        'F':{'TTC','TTT'},
        'Y':{'TAC','TAT'},
        '_':{'TAA'},
        'C':{'TGC'},
        'W':{'TGG'},   
	}


AA_table2 = {
            'I':{'ATC','ATT'},
            'M':{'ATG'},
            'T':{'ACC','ACG'} ,
        'N':{'AAC','AAT'},
        'K':{'AAA','AAG'},
        'S':{'AGC','AGT','TCT'},
        'R':{'CGC','CGT'},
        'L':{'CTG','CTT','TTA','TTG'},
        'P':{'CCA','CCG','CCT'},
        'H':{'CAC','CAT'},
        'Q':{'CAA','CAG'},
        'V':{'GTC','GTG','GTT'},
        'A':{'GCA','GCC','GCG'},
        'D':{'GAC','GAT'},
        'E':{'GAA','GAG'},
        'G':{'GGC','GGT'},
        'F':{'TTC','TTT'},
        'Y':{'TAC','TAT'},
        '_':{'TAA','TAG','TGA'},
        'C':{'TGC','TGT'},
        'W':{'TGG'},
        }



def getseq(file,line,column):
	x = open(file,'r')
	for y, z in enumerate(x):
		if y == line-1:
			string = z.split()[column-1]
	x.close()
	return string

def seq_from_gene(seq,table=codon_table):
#the table below contains the AA-codon pairs for the E.Coli K12 strain
#quick check indicates that table seem OK
  aa_sequence = ""
  if len(seq)%3 ==0:
  	for i in range(0,len(seq), 3):
  		codon = seq[i:i+3]
  		aa_sequence+=table[codon]
  else:
  	print('WARNING: genetic sequence not divisable by 3!')
  return aa_sequence
   
def pick_alternative_codon(seq, position_index,tableC=codon_table,tableA=AA_table,alternative=False):
  start_index = int(position_index)
  codon = seq[start_index:start_index+3]
  print('codon is:', codon)
  aa = deepcopy(tableC[codon])
  
  #if aa == "M" or aa == "C" or aa =="_" or aa =="W" or aa == "T" or aa == "K" or aa == "S" or aa == "L" or aa == "P" or aa == "V" or aa =="A":
  #	returnseq = deepcopy(seq)
  #	pass
  #else:

  possible_new_codons = deepcopy(tableA[aa])
  if alternative == True:
    possible_new_codons = deepcopy(AA_table2[aa])
  if codon in possible_new_codons and len(possible_new_codons) is not 1 :
    possible_new_codons.remove(codon)
  print('index is',start_index,'current codon is:',codon,'AA is:',aa,'choose from:',possible_new_codons)
  new_codon=rd.choice(tuple(possible_new_codons))
  print('Picked as a new codon:',new_codon)
  newseq = [c for c in seq]
  print('newseq check:', seq[start_index:start_index+3])
  newseq[start_index:start_index+3]=new_codon[:]
	#newseq = "".join(new_codon[i-start_index] if i == (start_index or start_index+1 or start_index +2) else c for i,c in enumerate(seq))
  returnseq = "".join(newseq)
  print('newseq doublecheck:',str(newseq[start_index:start_index+3]))
  return returnseq

def pick_start_numbers(base_range):
    start_list = np.array([])
    for i in base_range:
        if i % 3 == 0:
            start_list = np.append(start_list,(i-3))
        elif i % 3 == 2:
            start_list = np.append(start_list,(i-2))
        elif i % 3 == 1:
        	start_list = np.append(start_list,(i-1))
    print('start list is:',start_list)
    unique_list=np.unique(start_list)
    return unique_list

#testcase:
#aa_numbers_to_change = [ 15,24,25,28,36,40,44,52,54,55,89,94,100,117,119,123,125,127, \
#                        151,168,176,220,224,229,234,242,245,254,265,276,279,297,302,317, \
#                        331,345,350,373,388,392,403,409,415,422,437,444,451,452,458,463,\
#                        469,485,490,491,505,512,529,531,533,537,542,550,549,548,560,558,\
#                        565,563,575,581,593,610,627,634,638,649,647,669,685,701,709,737,\
#                        753,770,769,793,817]

#aa_numbers_to_change = [15,24,25,37,55,52,57,67,89,94,116,119,145,156,157,166,168,183,190,189,206,207,208,230,238,242,\
 #                       271,276,295,311,309,313,332,343,347,355,358,365,386,390,396,414,419,430,431,439,447,469,477,\
 #                       487,512,516,522,529,537,542,544,566,569,585,606,627,643,665,701,706,710,709,718,739,761,763,\
 #                       768,790,793,802]

#aa_numbers_to_change = [15,24,25,37,38,40,50,55,52,79,105,115,113,121,119,127,154,163,179,190,192,216,218,240,244,\
 #                       262,264,271,179,297,311,312,313,323,327,329,335,349,355,363,369,371,382,385,388,390,397,409,419,\
  #                      437,440,445,450,449,451,499,535,533,540,543,580,585,595,634,643,655,685,710,717,735,747,761,790]

#aa_numbers_to_change = [32,33,62,279,748,718]
#aa_numbers_to_change = [33,87,126,138,153,169,174,184,189,208,222,225,255,259,264,268,277,288,307,332,354,361,370, 382,395,414,439,791,777,767,748,739,735,706,719,718,665,673,687,693,617,635,640,648,569,598,541,455,494,497,499]

#aa_numbers_to_change = [168, 533]
#codon_numbers_to_change  = [3*aa_num for aa_num in aa_numbers_to_change]


codon_numbers_to_change = range(1744, 1744+30+1)
#codon_numbers_to_change2 = range(507,507+40+1)
#codon_numbers_to_change3  = range(627,627+28+1)
#codon_numbers_to_change4 = range(909,909+34+1)
#codon_numbers_to_change5 = range(1951,1951+41+1)

#codon_numbers_to_change.extend(codon_numbers_to_change2[:])
#codon_numbers_to_change.extend(codon_numbers_to_change3[:])
#codon_numbers_to_change.extend(codon_numbers_to_change4[:])
#codon_numbers_to_change.extend(codon_numbers_to_change5[:])

start_numbers = pick_start_numbers(codon_numbers_to_change)

#nupy_gene = getseq('tmp2.txt',1,1)
nupy_gene = getseq(path_file,1,1)
nupy_aa   = seq_from_gene(nupy_gene)

for number in start_numbers:
	nupy_gene = pick_alternative_codon(nupy_gene,number,alternative=False)

optimized_aa = seq_from_gene(nupy_gene)



if nupy_aa == optimized_aa:
	print('Sequence modification successful! Continuing.')

f = open('./optimized_genes/'+path_file[:-6]+"_CAI_fasta.txt","w+")
counter=0
for nucleotide in nupy_gene:
    f.write(nucleotide)
    counter+=1
    if counter % 60 == 0:
    	f.write('\n')
f.close()

f = open('./optimized_genes/'+path_file[:-6]+"_CAI.txt","w+")
counter=0
for nucleotide in nupy_gene:
    f.write(nucleotide)
f.close()



       
      #  'I':{'ATA','ATC','ATT'}, 
	  #  'M':{'ATG'},
	  #  'T':{'ACA','ACC','ACG','ACT'} ,
      #  'N':{'AAC','AAT'},
      #  'K':{'AAA','AAG'},
      #  'S':{'AGC','AGT','TCA','TCC','TCG','TCT'},
      #  'R':{'AGG','AGA','CGA','CGC','CGG','CGT'},
      #  'L':{'CTA','CTC','CTG','CTT'},
      #  'P':{'CCA','CCC','CCG','CCT'},
      #  'H':{'CAC','CAT'},
      #  'Q':{'CAA','CAG'},
      #  'V':{'GTA','GTC','GTG','GTT'},
      #  'A':{'GCA','GCC','GCG','GCT'},
      #  'D':{'GAC','GAT'},        
      #  'E':{'GAA','GAG'},
      #  'G':{'GGA','GGC','GGG','GGT'},
      #  'F':{'TTC','TTT'},
      #  'L':{'TTA','TTG'},
      #  'Y':{'TAC','TAT'},
      #  '_':{'TAA','TAG','TGA'},
      #  'C':{'TGC','TGT'},
      #  'W':{'TGG'},
            
