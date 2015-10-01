#!/usr/bin/env python

import os, sys
import scipy as SP
import h5py

#usage
def usage():
	print """ This script inputes missing values based on the mean of non-missing genotypes .
		Usage:
		
		python geno_preprocessing.py <chr1.hf5> """

#impute missing genotypes 
def impute_missing(matrix):

	'''this function imputes missing genotypes based on non-missing ones 
		for the same SNP across all the samples N '''

	for i in xrange(matrix.shape[1]):
		Imiss= SP.isnan(matrix[:,i])
		if not Imiss.any():
			continue
		else:
			mean = matrix[~Imiss,i].mean()
			matrix[Imiss,i] = mean
	return 'Done'

#calculate kinship 
def kinship(matrix, genofile):

	'''this function calculates matrix of correlation N x N, using
		matrix of N samples x G genotypes ''' 

	K = SP.dot(matrix,matrix.T)
	##populate Geno dataset with kinship matrix
	Kpopshape=K.shape
	kinship=genofile.create_dataset('genotype/Kpop',Kpopshape,dtype="float64")
	kinship[...]=K[:]
	return 'Done'

if __name__ == '__main__':
        
        if len(sys.argv[1:])!=1:
                sys.stderr.write("ERROR: geno file not provided\n")
                sys.stderr.write("Usage: geno_preprocessing.py hdf5_file\n")
                sys.exit(1)

	#check argument
	if len(sys.argv[1:]) != 1:
	        usage()
        	sys.exit(1)

	geno = sys.argv[1]

	if os.path.isfile(geno) != True:
        	sys.stderr.write("ERROR: file "+geno+" not found\n")
        	sys.exit(1)

	geno = h5py.File(geno)
	geno_matrix = geno['genotype/matrix'][:]
	impute_missing(geno_matrix)
	kinship(geno_matrix, geno)
	print 'check if keys are ok'
	for i in geno.keys():
    		for z in geno[i].keys():
        		print z

	geno.close()
