#!/usr/bin/env python

import sys, os
import peer
import h5py
import scipy as sp


def usage():
	print ''' 
This script runs PEER (https://github.com/PMBio/peer/wiki) on a matrix of normalised/not-normalised gene expression counts
with N samples and G genes.

Usage:
runpeer.py <pheno.filtered.hdf5> <hidden_k> <n_iterations> <peer_residuals.hd5> <covariates_sorted.hdf5>

NB:<covariates_sorted_hdf5> is optional. TODO: change make file.

'''


def runpeer(phenotype,K,iterations,cov=None):
	model=peer.PEER()
	model.setPhenoMean(phenotype)
	if cov is not None:
		model.setCovariates(cov)
	model.setNk(K)
	model.setNmax_iterations(iterations)
	model.update()
	residuals=model.getResiduals()
	return residuals

if __name__ == '__main__':

	cov_set=False

	if len(sys.argv[1:]) < 4:
		usage()
        	sys.stderr.write("\nERROR: missing parameter\n")
        	sys.exit(1)
	
	if len(sys.argv[1:]) == 4:
		print '\nRunning PEER without known covariates in the model.\n'
		pheno,hidden_k,n_iterations,outfile = sys.argv[1:]
	elif len(sys.argv[1:]) == 5:
		print '\nRunning PEER with known covariates in the model.\n'
		pheno,hidden_k,n_iterations,outfile,covariates = sys.argv[1:]
		if os.path.isfile(covariates) !=True:
			sys.stderr.write('\nERROR: file'+covariates+' not found\n')
			sys.exit(1)
		cov=h5py.File(covariates,'r')
		cov_matrix=cov['covariates'][:] #get matrix of covariates N x K. Catch warnings here
		cov_set=True

	if os.path.isfile(pheno) != True:
		sys.stderr.write("\nERROR: file "+pheno+" not found\n")
		sys.exit(1)

	pheno=h5py.File(pheno, 'r')

	#run peer
	sample_ID = pheno['phenotype/row_header/sample_ID'][:] #get samples IDs; catch warnings here
	gene_ID = pheno['phenotype/col_header/phenotype_ID'][:] #get genes IDs; catch warnings here
	phenotype = pheno['phenotype/Ytransformed'][:] #get matrix of gene counts; catch warnings here
	samples = phenotype.shape[0]
	#exit if the number of hidden confounding factor selected is >=25% of samples used in the analysis!
	threshold = (25.0*samples)/100
        #select hidden confounding 
	hidden_k = int(hidden_k)

	if hidden_k > threshold:
		sys.stderr.write('\nERROR: Number of hidden factors chosen is above 25% of the number of samples. Please select a lower value.\n')
		sys.exit(1)
			 
	#iterations and outfile
        n_iterations=int(n_iterations)
        outfile=h5py.File(outfile,'w')
	#apply PEER model
	if cov_set==False:
		residuals = runpeer(phenotype,hidden_k,n_iterations)
	else:
		residuals = runpeer(phenotype,hidden_k,n_iterations,cov=cov_matrix)
	#write within out file
	dset=outfile.create_dataset('phenotype',residuals[:].shape, dtype='float64')
	dset[...]=residuals[:]
	dset2=outfile.create_dataset('row_header/sample_ID',sample_ID.shape,dtype='S1000')
	dset2[...]=sample_ID
	dset3=outfile.create_dataset('col_header/phenotype_ID',gene_ID.shape,dtype='S1000')
	dset3[...]=gene_ID

	outfile.close()

	sys.exit(0)
