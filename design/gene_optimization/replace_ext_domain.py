import numpy
import sys
import os

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
seqname = seqfilename.split('.')
outputfilename = seqname[0] + '_replaced_ext_domain.fasta'

sequence_name = getseq(seqfilename,1,1)

os.system("grep -v '>' %s > tmp.txt" % (seqfilename))
os.system("tr -d '\n' < tmp.txt > tmp2.txt")

sequence = getseq('tmp2.txt',1,1)

seqname = seqfilename.split('.')
#print(seqname)






ext_domain_seq = 'ACCTGAAATCTTCAATTGGTAAAGACTTCACCTTAGCGCAGGTGAACAAAACGACCAACAACGCAAGCATCATTGAGGCCATTCCTTATTATAATCTGGTCAAGAAGTCTAGTAATAAACGCCGTAGTTATAACGAGGATGACCCGCGCTTAGGCGAGAACGCGTCGTCTGCCAAAGCCAAAAAATTTCAGCTGGATCGCGAACAGTTAAATTCAGCGGGTAAAGATGACTTTGACATCCAGTTATTATCGCCGCACGAGTCGGTGTTCATGATCATGCTGAACGCCGGTAACTTAGATGCGCGTTTAGACGCAAATAATAAAGTCGAAAGCTTAAATTCAAACGATGCGCAGCGCAGCAACGGCGAAAACCCGAAACCGCGCGCGATGACTCCGCTGACCATCAGCAGCGAAAAAGATTCTATCAGCTCAAGCATGAAAATCGATCCGGCGAAACATCTGACTGAGCCGAGCCAGGACATTGAACAGAAAAAACTGACCACCACCTGGAAGTCATCAGATGAAAAGCCGGTTAGTGTGCCGAAATCGACCAACTCATCGGAAGAATGCTAA'
ext_domain_list = list(ext_domain_seq)

end_seq = len(sequence)-1
counter = 0
for i in reversed(range(len(ext_domain_list))):
	sequence[end_seq-counter]=ext_domain_list[i]
	counter+=1


with open(outputfilename,"w") as file:
    file.write("> %s\n" % (seqname[0]))
    counter = 0.
    for cod in sequence:
        file.write("%s" % cod)
        counter+=1
        if counter % 60 == 0:
            file.write("\n")


