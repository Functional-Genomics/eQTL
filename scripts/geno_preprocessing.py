#!/usr/bin/env python

""" This script inpute missing values by mean of non-missing ones.
	usage> python geno_preprocessing.py chr1.hf5 """

import os, sys
limix_path=''
sys.path.insert(0,limix_path)
import scipy as SP
import h5py

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

	geno = sys.argv[1]

	if os.path.isfile(geno)!=True:
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
