import sys
import numpy as np

non_def_pairs = str(sys.argv[1])

if non_def_pairs == True:
	sigma_1 = sys.argv[2]
	sigma_2 = sys.argv[3]
	entity_1 = sys.argv[4]
	entity_2 = sys.argv[5]
	sigma = (sigma_1+sigma_2)/2
else:
	sigma = 0.6 #default for amino acids

'''
Built originally by A. Jansen, MSc, further edits by H.W. de Vries (2020).
This code creates non-bonded interactions for use in GROMACS. These are stored in 
a table file and comprise electrostatics and LJ-type hydrophobic potentials.
See also: http://manual.gromacs.org/documentation/current/reference-manual/special/tabulated-interaction-functions.html
http://www.gromacs.org/Documentation/How-tos/Tabulated_Potentials#Non-bonded_tabulated_interactions
These potentials follow those outlined in the original 1-BPA paper, see
Ghavami, et al. BPJ for details on the model.

Note that for electrostatics, a monovalent salt concentration of 150mM
is used as a default (see ionCalc) together with a temp of 300K.

Inputs:
- non_def_pairs = boolean(True or False): if set to False, no further inputs
                  are required, and the code assumes that amino acids are used.
         	  This is useful for generating tables with different ionic strength

- sigma_1 (opt) = float(in nm), sets the size (sigma) of particle one.
- sigma_2 (opt) = as above, for particle two.
- entity_1 (opt) = string (pref. full-caps). Name of the type of bead associated
   	 	   with sigma_1. This should match the name of an energy group
                   involved in the interaction for which the table is made.
- entity_2 (opt) = string (pref. full-caps). As above. 

Outputs:
- table.xvg  = GROAMCS-readable file that is multiplied with the parameters
               in forcefield.itp. Note that vdw-type=user and coul-type=user
               should be set in the .mdp file!
- table_<entity_1>_>entity2>.xvg (opt)= GROMACS readable file that is multiplied
               with the parameters in forcefield.itp. Same conditions apply,
               but in this case the two entities need their own index groups
               AND need to be mentioned as specific table-pairs under 
               energygrp_table, as in: <entity1> <entity2>. 
'''

def ioncalc(c=150,T=300):
	'''
	Calculates the kappa factor in the Debye-Huckel term in the electrostatics.
	Assumes a standard value of 150mM ionic strength (i.e. 150mM monovalent salt) and T=300K
	'''

	e_0 = 8.85418781*10**-12							#Calculate kappa as a function of concentration and T
	k_b = 1.38064853*10**-23
	e   = 1.60217662*10**-19
	N_A = 6.02214076*10**23
	e_r = 78

	I   = c
	d   = ((e_0*e_r*k_b*T)/(2*N_A*(e**2)*I))**(0.5)
	k   = (d*10**9)**-1
	
	return k

def f(r,k):
	S_s   = 1.0

	z     = 0.25d
	frac1 = (r**2)/(z**2)								#Calculate sigmoidal function
	frac2 = (np.exp(r/z)) / ((np.exp(r/z)-1)**2)
	e_r   = S_s*(1 - frac1*frac2)

	f     = (e_r*r)**(-1) * np.exp(-k*r)				#Calculate modified Coulomb potential

	return f

def df(flist):
	dy    = np.diff(flist)/0.002						#Calculate derivative of f using numpy.
	dy    = np.ndarray.tolist(dy)

	dy.append(dy[len(dy)-1])
	dy[149] = 0											#Repairs a mistake made by np.diff algorithm.

	dz = []
	for i in dy:
		dz.append(-i)

	return dz

# FUNCTIONS FOR ALI-8-6

def g_ali(r,sigma,g_trans,dg_trans):
	if r < 0.3:
		return (0.3*dg_trans+g_trans)-dg_trans*r 		#Linear function with slope of g(0.3)
	elif r < sigma:
		return -(4.0/3.0)*(sigma/r)**6 + (1.0/3.0)
	else:
		return -(sigma/r)**8

def dg_ali(r,sigma,dg_trans):
	if r < 0.3:
		return dg_trans									# = Slope of g(0.3) = dg(0.3)
	elif r < sigma:
		return -(8*sigma**6)/(r**7)
	else:
		return -(8*sigma**8)/(r**9)

def h_ali(r,sigma,h_trans,dh_trans):
	if r < 0.3:
		return (0.3*dh_trans+h_trans)-dh_trans*r 		#Linear function with slope of h(0.3)
	else:
		return (sigma/r)**8

def dh_ali(r,sigma,dh_trans):
	if r < 0.3:
		return dh_trans									# = Slope of h(0.3) = dh(0.3)
	else:
		return (8*sigma**8)/(r**9)

###################################################################
# CREATE table.xvg
###################################################################

print('Building table files...')

k = ioncalc()

#GENERATE rlist

rlist = []
for i in xrange(0,3001):
	rlist.append(0.002*i)

#CALCULATE ELECTROSTATICS: f(r), df(r)

flist = []

for r in rlist:
	if r < 0.3:
		flist.append(0)									#We set f(r) to zero for r < 0.3.
	else:
		flist.append(f(r,k))

dflist = df(flist)

#CALCULATE HYDROPHOBICS: g_ali(r), dg_ali(r)

glist = []; dglist = []

g_trans  = g_ali(0.3,sigma,1,1)							#For calculating our line for r < 0.3.
dg_trans = dg_ali(0.3,sigma,1)

for r in rlist:
	glist.append(g_ali(r,sigma,g_trans,dg_trans))
	dglist.append(dg_ali(r,sigma,dg_trans))

#CALCULATE HYDROPHOBICS: h_ali(r), dh_ali(r)

hlist = []; dhlist = []

h_trans  = h_ali(0.3,sigma,1,1)							#For calculating our line for r < 0.3.
dh_trans = dh_ali(0.3,sigma,1)

for r in rlist:
	hlist.append(h_ali(r,sigma,h_trans,dh_trans))
	dhlist.append(dh_ali(r,sigma,dh_trans))

#PRINTING TABLE.XVG

if non_def_pairs==True:
	ff = open("table_%s_%s.xvg" % (entity_1,entity_2),"w+")
else:
	ff = open("table.xvg","w+")

for i in xrange(0,len(rlist)):
     ff.write("%E  %E  %E  %E  %E  %E  %E \n" % (rlist[i],flist[i],dflist[i],glist[i],dglist[i],hlist[i],dhlist[i]))

ff.close()
