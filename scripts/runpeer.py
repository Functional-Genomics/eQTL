#!/usr/bin/env python

import sys, os
#sys.path.append('./../')
# this is really necessary?
#peer_path='/nfs/research/stegle/software/opt/peer/build/python' #add this to the install.sh script?
#sys.path.insert(0,peer_path)
import peer
import h5py
import scipy as sp


def usage():
	print ''' 
This script runs PEER (https://github.com/PMBio/peer/wiki) on a matrix of normalised/not-normalised gene expression counts
with N samples and G genes.

Usage:
runpeer.py <pheno.filtered.hdf5> <hidden_k> <n_iterations> <peer_residuals.hd5>

TODO: implement covariates in the model. It's one line of code but I need to test it first. '''


def runpeer(phenotype,K,iterations):
	model=peer.PEER()
	model.setPhenoMean(phenotype)
	model.setNk(K)
	model.setNmax_iterations(iterations)
	model.update()
	residuals=model.getResiduals()
	return residuals

if __name__ == '__main__':

	if len(sys.argv[1:]) < 4:
        	sys.stderr.write("ERROR: missing parameter\n")
        	usage()
        	sys.exit(1)

	pheno,hidden_k,n_iterations,outfile = sys.argv[1:]
	
	if os.path.isfile(pheno) != True:
		sys.stderr.write("ERROR: file "+pheno+" not found\n")
		sys.exit(1)

	pheno=h5py.File(pheno, 'r')
	hidden_k = int(hidden_k)
	n_iterations=int(n_iterations)
	outfile=h5py.File(outfile,'w')
	
	#run peer
	phenotype = pheno['phenotype/Ytransformed'][:] #catch warnings here
	residuals = runpeer(phenotype,hidden_k,n_iterations)
	#write within out file
	dset=outfile.create_dataset('phenotype',residuals[:].shape, dtype='float64')
	dset[...]=residuals[:]

	outfile.close()

	sys.exit(0)
