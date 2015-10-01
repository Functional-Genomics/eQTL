#!/usr/bin/env python

import os, sys
import scipy as SP
import h5py

#usage
def usage():
	print """ This script inputes missing values based on the mean of non-missing genotypes .
Usage: geno_preprocessing.py hdf5_file"""

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
def corr_matrix(matrix, genofile):

	'''this function standardises the geno data by centering SNPs to mean 0, with unit variance. 
	Then it calculates matrix of correlation N x N, using matrix of N samples x G genotypes ''' 
	matrix -= matrix.mean(0)
	matrix /= matrix.std(0)
	K = SP.dot(matrix,matrix.T)
	##populate Geno dataset with kinship matrix
	Kshape=K.shape
	kpop=genofile.create_dataset('genotype/Kpop',Kshape,dtype="float64")
	kpop[...]=K[:]
	return 'Done'

if __name__ == '__main__':
        

	#check argument
	if len(sys.argv[1:]) != 1:
                sys.stderr.write("ERROR: geno hdf5 file not provided\n")
	        usage()
        	sys.exit(1)

	geno = sys.argv[1]

	if os.path.isfile(geno) != True:
        	sys.stderr.write("ERROR: file "+geno+" not found\n")
        	sys.exit(1)

	geno = h5py.File(geno)
	geno_matrix = geno['genotype/matrix'][:]
	impute_missing(geno_matrix)
	corr_matrix(geno_matrix, geno)
	print 'check if keys are ok'
	for i in geno.keys():
    		for z in geno[i].keys():
        		print z

	geno.close()
