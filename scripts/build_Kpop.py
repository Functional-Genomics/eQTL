#!/usr/bin/env python

import os, sys
import numpy as np
import h5py

def usage():

        print ''' This script builds the final genetic kinship matrix by summing up all the chromosomal kinship matrices given.

Usage: build_Kpop.py <Kpop.hdf5> <chr1.hdf5> <chr2.hdf5> ... <chr22.hdf5> '''

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
	return Kpop,samples

if __name__ == "__main__":
	#check arguments
	if 'Kpop' not in sys.argv[1]:
		sys.stderr.write("ERROR: Kpop out file not provided\n")
		usage()
		sys.exit(1)
	elif len(sys.argv[1:]) < 2 :
		sys.stderr.write("ERROR: chromosome not provided\n")
		sys.exit(1)

	#arguments
	Kpopfile = h5py.File(sys.argv[1],'w')
	chr = sys.argv[2:]

	#check files exist
	for file in chr:
		if os.path.isfile(file) != True:
			Kpopfile.close()
			sys.stderr.write("ERROR: file "+file+" not found\n")
			sys.exit(1)

 	#populating Kpop matrix
	samples = ''
	Kpop = ''
	for file in chr:
    		X=h5py.File(file,'r')
		matrix = X['genotype/Kpop'][:]
		Kpop,samples = build_kpop(matrix,Kpop,samples)
		X.close()

	#output Kpop in a file
	kinship=Kpopfile.create_dataset('Kpop',Kpop.shape,dtype="float64")
	kinship[...]=Kpop[:]
	Kpopfile.close()

	sys.exit(0)

