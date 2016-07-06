#!/usr/bin/env python

import os, sys
import scipy as SP
import pandas as pd
import h5py

#usage
def usage():
	print """ This script inputes missing values based on the mean of non-missing genotypes and calculates genetic kinship. If skip_kinship parameter is provided it generates a Kpop matrix of NaN.

Usage: geno_preprocessing.py <hdf5_file> 

Optional parameters: <skip_kinship> """

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
	return matrix

#calculate kinship 
def corr_matrix(matrix,skip_kinship):

	'''if calculate_kinship is False, this function standardises the geno data by centering SNPs to mean 0, with unit variance. Then it calculates matrix of correlation N x N, using matrix of N samples x G genotypes. If skip_kinship is True, it generates a N x N matrix of NaN '''
	if skip_kinship == False: 
		matrix -= matrix.mean(0)
		matrix /= matrix.std(0)
		nan=SP.isnan(matrix)
		matrix[nan]=0.0
		K = SP.dot(matrix,matrix.T)
	else:
		#if kinship is not requested
		shape=matrix.shape[0]
		K = SP.empty((shape,shape))
		K[:] = SP.NAN
	return matrix,K


if __name__ == '__main__':
        
	skip_kinship=False

	#check argument
	if len(sys.argv[1:]) < 1:
		usage()
                sys.stderr.write("ERROR: missing parameter\n")
        	sys.exit(1)

	geno = sys.argv[1]

	if os.path.isfile(geno) != True:
        	sys.stderr.write("ERROR: file "+geno+" not found\n")
        	sys.exit(1)
	
	#import genotype
	genofile = h5py.File(geno)
	geno_matrix = genofile['genotype/matrix'][:]
	#substitute missing values with mean of not-missing 
	genofile['genotype/matrix'][:] = impute_missing(geno_matrix)
	#standardise matrix of SNPs and get Kpop or generate a NaN Kpop matrix if requested
	if len(sys.argv[1:]) == 2:
		print 'generating a matrix of NaN as Kpop'
		skip_kinship = True
		matrix,K = corr_matrix(geno_matrix,skip_kinship)
	else:
		matrix,K = corr_matrix(geno_matrix,skip_kinship)
	#substitute matrix of imputed values with standardise matrix
	genofile['genotype/matrix'][:]=matrix[:]
	#set a dataset with var_names
	chrom = pd.Series(genofile['genotype/col_header/chrom'][:]).astype(int).astype(str)
	pos = pd.Series(genofile['genotype/col_header/pos'][:]).astype(str)
	var_names = chrom.str.cat(pos,sep='_').values
	#populate dataset with Kpop
        Kshape=K.shape
        kpop=genofile.create_dataset('genotype/Kpop',Kshape,dtype="float64")
        kpop[...]=K[:]
	dset = genofile.create_dataset('genotype/col_header/var_names',data=var_names.tolist())
	print 'check if file {0} has 4 keys: Kpop, col_header, row_header, matrix\n'.format(geno)
	for dset in genofile.keys():
    		for key in genofile[dset].keys():
        		print key
	print 'OK...all the keys found'
	
	
	genofile.close()

	sys.exit(0)
