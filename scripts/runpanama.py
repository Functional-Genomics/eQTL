#!/usr/bin/env python

import sys, os
import scipy as sp
import h5py
import os
# import LIMIX
# here I'm afraid that the name of the modules changed in the meanwhile. Check the latest version of limix TODO: check
import limix.modules.panama as panama 

def usage():
        print """
This module runs  the PANAMA model available within LIMIX, to calculate matrix of relatedness by using genetic kinship and the phenotype matrix. User can define the number of ranks (K) to use. Selection of most variable genes in based on signal-to-noise ratio calculation for each gene.

TODO: Add other methods to select genes from pheno matrix.

Usage:
runpanama.py <Kpop.hdf5> <pheno.filtered.hf5> <hidden_k> <snr_threshold> <fileout.hdf5> """ 

def standardise_kinship(matrix):
	kinship /= sp.diag(matrix).mean()
	kinship += 1e-4 * sp.eye(kinship.shape[0])
	return kinship

def select_genes_on_SNR(phenotype, t):
	#calculate the signal-to-noise ratio for each gene of the matrix
	SNR = (phenotype.mean(0)/phenotype.std(0))
	#get the percentile of the SNR distribution based on threshold user defined
	perc = sp.percentile(sp.sort(SNR),t)
	#take indexes of those genes with SNR >=t
	index = sp.where(SNR >= perc)[0]
	#take genes from phenotype matrix with SNR >= t
	Y = sp.take(phenotype,index, axis=1)
	return Y

def panama(phenotype, Kpop, hidden_k):
	#build cumulative distribution of variance explained by the selected genes
	cum_var = panama.PC_varExplained(phenotype)
	#train panama with user defined ranks
	p = panama.PANAMA(Y=phenotype, Kpop = Kpop)	
	p.train(rank=hidden_k)
	# retrieve stuff
	Kpanama = p.get_Kpanama()
	Ktot    = p.get_Ktot()
	vComp   = p.get_varianceComps()
	return Kpanama,Ktot,vComp

if __name__ == '__main__':

	sp.random.seed(0) #TODO: not sure why they got this in the code (is to set the same starting seed). Maybe for panama training. I have to check this.
	if len(sys.argv[1:]) < 5:
        	sys.stderr.write('ERROR: missing parameters\n')
        	usage()
        	sys.exit(1)

	Kpop,pheno,hidden_k,snr_threshold,fileout = sys.argv[1:]

	if os.path.isfile(Kpop) != True:
	        sys.stderr.write('ERROR: '+Kpop+' not found\n')
	        sys.exit(1)
	elif os.path.isfile(pheno) != True:
	        sys.stderr.write('ERROR: '+pheno+' not found\n')
	        sys.exit(1)

	#arguments
	Kpop = h5py.File(Kpop, 'r')
	pheno = h5py.File(pheno, 'r')
	hidden_k = int(hidden_k)
	snr_threshold = int(snr_threshold)
	fileout = h5py.File(fileout, 'w')

	#read matrix
	phenotype = pheno['phenotype/Ytransformed'][:] #catch warnings here
	kinship = Kpop['Kpop'][:] #catch warnings here	
	
	#apply func
	stdkinship = standardise_kinship(kinship)
	selected_phenotypes = select_genes_on_SNR(phenotype, snr_threshold)
	print "{0} genes selected, featuring a signal-to-noise ratio >= {1} percentile".format(selected_phenotypes.shape[1],snr_threshold)
	Kpanama,Ktot,vComp = panama(selected_phenotypes[:], stdkinship[:], hidden_k)
	print "Variance composition:\n{0}".format(vComp)
	#writing into fileout	
	Ktot_shape=Ktot[:].shape
	dset=filout.create_dataset('Ktot', Ktot_shape, dtype='float64')
	dset[...]=Ktot[:]

	fileout.close()
	print "all done."
	sys.exit(0)
