#!/usr/bin/env python

import os, sys
import numpy as np
import h5py

def usage():

        print ''' This script builds the final genetic kinship matrix (Kpop.hdf5) by summing up all the chromosomal kinship matrices given.

Usage: build_Kpop.py <Kpop.hdf5> <chr1.hdf5> [<chr2.hdf5> ... ]'''

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
        if len(sys.argv[1:]) == 0 :
		sys.stderr.write("ERROR: missing parameters\n")
                usage()
		sys.exit(1)
	elif len(sys.argv[1:]) < 2 :
		sys.stderr.write("ERROR: full set of parameters not provided \n")
		usage()
		sys.exit(1)

	#arguments
	chr = sys.argv[2:]

	#check if files exist first
	for file in chr:
		if os.path.isfile(file) != True:
			sys.stderr.write("ERROR: file "+file+" not found\n")
			sys.exit(1)
        
	#open Kpop out file                
	Kpopfile = h5py.File(sys.argv[1],'w')
 	#populating Kpop matrix
	samples = ''
	Kpop = ''
	for file in chr:
    		X=h5py.File(file,'r' ) #catch warning if file is corrupted
		matrix = X['genotype/Kpop'][:]
		Kpop,samples = build_kpop(matrix,Kpop,samples)
		X.close()

	#output Kpop in a file
	kinship=Kpopfile.create_dataset('Kpop',Kpop.shape,dtype="float64")
	kinship[...]=Kpop[:]
	Kpopfile.close()

	sys.exit(0)

