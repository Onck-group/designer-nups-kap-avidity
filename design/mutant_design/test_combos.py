import numpy
import math


#AA =             ['A','R','N' , 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
#hydrophobicity = [0.7,0.0,0.33,0.35,0.68,0.41,0.35,0.64,0.53,0.98,1.00,0.15,0.78,1.00,0.75,0.45,0.51,0.96,0.82,0.94]

#AA =             ['A','R','N' , 'D', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
#hydrophobicity = [0.7,0.0,0.33,0.35,0.41,0.35,0.64,0.53,0.98,1.00,0.15,0.78,1.00,0.75,0.45,0.51,0.96,0.82,0.94]

AA =            ['A','R','N' , 'D'  , 'C', 'Q', 'E'  , 'G', 'H', 'I', 'L', 'K'  , 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
hydrophobicity= [0.7,0.0,0.33,0.0005,0.68,0.675,0.0005,0.41,0.53,0.98,1.00,0.0005,0.78,1.00,0.65,0.45,0.51,0.96,0.82,0.94]


for residue,scale in zip(AA,hydrophobicity):
    for residue2,scale2 in zip(AA,hydrophobicity):
        combined_value = scale+scale2
        if abs( combined_value - 1.41 ) < 0.05:
            print('MATCH! between:',residue,' and ',residue2, ', value is: ', combined_value)
