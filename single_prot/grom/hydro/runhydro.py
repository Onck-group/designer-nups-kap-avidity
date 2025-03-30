#!/bin/python

import os
import numpy

# PARAMETERS
##########################################################################################

# Temperature (C). We use T = 293 K as Ali also used this.
# The 300 vs 293 difference does not influence results in any way (I checked).
T      = 20				

# Average bead mass (Da) and average residue radius (A).
# See p. 25 Ghavami thesis and Zhang2002 where this value is coming from.
m      = 124
Br     = 5.1

# Number of frames to analyze.
frames = 500

# CLEANUP
##########################################################################################

os.system('rm -r frames      &> /dev/null')
os.system('rm -r frames_new  &> /dev/null')
os.system('rm -r outputs     &> /dev/null')

os.system('rm framelist.out  &> /dev/null')

os.system('mkdir frames')
os.system('mkdir frames_new')
os.system('mkdir outputs')

# REPAIR PARSING ERRORS IN MD_chain_hydro.pdb
##########################################################################################

def list2strlist(file):
	with open(file) as f:
		x = f.read().splitlines()
	return x

chain = list2strlist('./MD_chain_hydro.pdb')

f = open('MD_chain_hydro2.pdb','w+')

for line in chain:
	if len(line) < 78:
		f.write('%s\n' % (line))

	if len(line) >= 78:
		a = []
		for i in xrange(0,len(line)):
			if line[i] == ".":
				a.append(i)

		a1 = line[0:(a[0]+4)]
		a2 = line[(a[0]+4):(a[1]+4)]
		a3 = line[(a[1]+4):len(line)]

		f.write('%s %s %s\n' % (a1,a2,a3))

f.close()

# SPLIT .PDB FILE INTO INDIVIDUAL FRAMES
##########################################################################################

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

beads = file_len('chain_org2.pdb') - 6

os.system('split -l %s MD_chain_hydro2.pdb ./frames/frame -d' % (beads+6))

# CREATE FRAMELIST FOR FOR-LOOP BY PICKING FRAMES AT RANDOM
##########################################################################################

os.system('cd frames; ls > ./../framelist.out; cd ..')

framelist = list2strlist('framelist.out')

framelist = numpy.random.choice(framelist, frames)

# CONVERT FRAMES TO CORRECT FORMAT FOR HYDRO
##########################################################################################

def getcolumn(file,column,start,stop):
	q = []
	x = open(file,'r')
	for y, z in enumerate(x):
		if y >= start-1 and y <= stop-1:
			string = z.split()[column-1]
			q.append(float(string))
	x.close()
	return q

for i in framelist:
	
	x = getcolumn('./frames/%s' % (i), 6, 5, beads+4)
	y = getcolumn('./frames/%s' % (i), 7, 5, beads+4)
	z = getcolumn('./frames/%s' % (i), 8, 5, beads+4)

	f = open('%s' % (i),'w+')

	f.write('1.E-08,    !Unit of length for coordinates and radii, cm (10 A)\n')
	f.write('%s,        !Number of beads\n' % (beads))	

	for j in xrange(0,len(x)):
		f.write('%s \t %s \t %s \t %s\n' % (x[j], y[j], z[j], Br))

	f.close()

	os.system('mv %s frames_new' % (i))

# PRINT PARAMETER FILE FOR HYDRO++10
##########################################################################################

f = open('params.txt','w+')

f.write('whatever                 Title\n')
f.write('output                  filename for output files\n')
f.write('test.dat              Structural (bead coords) file\n')
f.write('12                              ICASE\n')
f.write('%s.                             Temperature, centigrade\n' % (T))
f.write('0.010                           Solvent viscosity\n')
f.write('%s.                          Molecular weight\n' % (beads*m))
f.write('0.760                           Specific volume of macromolecule\n')
f.write('1.0                             Solution density\n')
f.write('26,                 Number of values of H\n')
f.write('1.5e+7,             HMAX\n')
f.write('30,                 Number of intervals for the distance distribution\n')
f.write('80.E-8              RMAX\n')
f.write('10000,             (ONLY IF ISCA IS NOT ZERO) NTRIALS\n')
f.write('1                   IDIF=1 (yes) for full diffusion tensors\n')
f.write('*           End of file\n')

f.close()

# MAIN LOOP
##########################################################################################

for i in framelist:

	# GET FRAME

	os.system('cd frames_new; mv %s ./..; cd ..; mv %s test.dat' % (i,i))

	# RUN HYDRO++10 ON FRAME

	os.system('echo params.txt | ./hydro++10-lnx.exe')

	# CLEANUP

	os.system('mv output-res.txt ./outputs/output-res_%s.txt' % (i))
	os.system('rm params-* test.dat output*')

# ANALYSIS
##########################################################################################

def getfloat(file,line,column):
	x = open(file,'r')
	for y, z in enumerate(x):
		if y == line-1:
			number = z.split()[column-1]
	x.close()
	return float(number)

def stdev(A):
	N = len(A)
	mean = sum(A)/float(N)

	y = 0
	for i in A:
		y += (i-mean)**2

	return (y/float(N-1))**0.5

Rh_list = []

for i in framelist:
	val = getfloat('./outputs/output-res_%s.txt' % (i), 150, 2)
	
	if val != -1.0:
		Rh_list.append(val)

Rh_mean  = sum(Rh_list)/len(Rh_list)
Rh_stdev = stdev(Rh_list)

# WRITE HYDRO.OUT
##########################################################################################

f = open('hydro.out','w+')

f.write('%s        %s      nm\n' % (Rh_mean*10**7, Rh_stdev*10**7))
f.write('%s        %s       A\n' % (Rh_mean*10**8, Rh_stdev*10**8))
f.write('%s        %s      cm\n\n' % (Rh_mean, Rh_stdev))
f.write('rStokes(mean)        rStokes(stdev)\n\n')

for i in xrange(0, len(Rh_list)):
	f.write("%s    %.3E\n" % (framelist[i], Rh_list[i]))

f.close()

# CLEANUP
##########################################################################################

os.system('rm MD_chain_hydro.pdb   &> /dev/null')
os.system('rm MD_chain_hydro2.pdb  &> /dev/null')

os.system('rm params.txt           &> /dev/null')
os.system('rm framelist.out        &> /dev/null')

os.system('rm -r frames            &> /dev/null')
os.system('rm -r frames_new        &> /dev/null')

os.system('rm -rf outputs	   &> /dev/null')
