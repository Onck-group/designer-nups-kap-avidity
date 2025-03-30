import numpy
import sys


def getseq(file,line,column):
	'''
	reads in the sequence from a source file, assuming that everything is placed on one very long line.
	This might also work for FASTA etc, but requires a different way of being called. Use with caution!
	Credits to Anton Jansen for putting this black box together.
	'''
	x = open(file,'r')
	for y, z in enumerate(x):
		if y == line-1:
			string = z.split()[column-1]
	x.close()
	list_string=[cod for cod in string]
	return list_string


seqfilename = sys.argv[1]
sequence = getseq(seqfilename,1,1)
seqname = seqfilename.split('.')
#print(seqname)
outputfilename = seqname[0] + '_added_tags.fasta'

scissase_seq = 'ATGGGCCACCATCACCATCACCATCACCATGATTACGATATTCCAACGACCCTGGAAGTTCTGTTCCAGGGGCCC'
cys_term = 'TGCTAA'

sequence = list(scissase_seq)+sequence+list(cys_term)
#print(sequence)
with open(outputfilename,"w") as file:
    file.write("> %s\n" % (seqname[0]))
    counter = 0.
    for cod in sequence:
        file.write("%s" % cod)
        counter+=1
        if counter % 60 == 0:
            file.write("\n")

