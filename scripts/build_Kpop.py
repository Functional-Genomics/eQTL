#!/usr/bin/env python

import os, sys
import numpy as np
import h5py

def usage():

        print ''' This script builds the final genetic kinship matrix (Kpop.hdf5) by summing up all the chromosomal kinship matrices given.

Usage: build_Kpop.py <Kpop.hdf5> <samples.hdf5> <chr1.hdf5> [<chr2.hdf5> ... ] '''


def build_kpop(chrmatrix,Kpop,samples):
	if samples == '':
		samples = chrmatrix.shape
	else:
		pass
	if Kpop == '':
		Kpop=np.zeros(samples,dtype='float64')
		Kpop += chrmatrix
	else:
		Kpop += chrmatrix
	#normalise Kpop
	Kpop /= Kpop.diagonal().mean()
	return Kpop,samples

if __name__ == "__main__":
	#check arguments
        if len(sys.argv[1:]) == 0 :
		usage()
		sys.stderr.write("\nERROR: missing parameters\n")
		sys.exit(1)
	elif len(sys.argv[1:]) < 3 :
		usage()
		sys.stderr.write("\nERROR: full set of parameters not provided\n")
		sys.exit(1)

	#arguments
	chr = sys.argv[3:]

	#check if files exist first
	for file in chr:
		if os.path.isfile(file) != True:
			sys.stderr.write("\nERROR: file "+file+" not found\n")
			sys.exit(1)

        print "Opening input files..."
	#open Kpop out file                
	Kpopfile = h5py.File(sys.argv[1],'w')
	#open samples file
	hdf5_samples = h5py.File(sys.argv[2],'w')

 	#populating Kpop matrix and vector of samples
	samples = ''
	Kpop = ''
	samples_vector = 0
	for file in chr:
                print "."                                
    		X=h5py.File(file,'r' ) #catch warning if file is corrupted
		matrix = X['genotype/Kpop'][:]
		Kpop,samples = build_kpop(matrix,Kpop,samples)
		if samples_vector == 0:
			samples_vector = X['genotype/row_header/sample_ID'][:] #get samples from chr file
		X.close()

        #print "Creating kinship matrix..."
	print "\nCreating output file...\n"                           
	#output Kpop in hdf5 file
	kinship=Kpopfile.create_dataset('Kpop',Kpop.shape,dtype="float64")
	kinship[...]=Kpop[:]
	#output samples list into hdf5 file and also within Kpop file
	s=hdf5_samples.create_dataset('sample_ID', samples_vector.shape, dtype='S1000')
	s[...]=samples_vector[:]
	s1=Kpopfile.create_dataset('row_header/sample_ID',samples_vector.shape,dtype='S1000')
	s1[...]=samples_vector[:]
	s2=Kpopfile.create_dataset('col_header/sample_ID',samples_vector.shape,dtype='S1000')
	s2[...]=samples_vector[:]
	Kpopfile.close()
        hdf5_samples.close()
	print "\nCreating output file...done.\n"

	sys.exit(0)

