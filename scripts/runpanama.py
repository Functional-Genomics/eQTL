#!/usr/bin/env python

import sys
import scipy as sp
import h5py
import os
# import LIMIX
# here I'm afraid that the name of the modules changed in the meanwhile. Check the latest version of limix TODO: check
import limix.modules.panama as panama 
from limix.stats.pca import *

def usage():
        print """
This module runs  the PANAMA model available within LIMIX, to calculate matrix of relatedness by using genetic kinship and the phenotype matrix. User can define the number of ranks (K) to use. Selection of most variable genes in based on signal-to-noise ratio calculation for each gene.

TODO: Add other methods to select genes from pheno matrix.

Usage:
runpanama.py <Kpop.hdf5> <pheno.filtered.hf5> <hidden_k> <snr_threshold> <fileout.hdf5> """ 

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

if __name__ == '__main__':

	sp.random.seed(0) #TODO: not sure why they got this in the code (is to set the same starting seed). Maybe for panama training. I have to check this.
	if len(sys.argv[1:]) < 5:
		usage()
        	sys.stderr.write('\nERROR: missing parameters\n')
        	sys.exit(1)

	Kpop,pheno,hidden_k,snr_threshold,fileout = sys.argv[1:]

	if os.path.isfile(Kpop) != True:
	        sys.stderr.write('\nERROR: '+Kpop+' not found\n')
	        sys.exit(1)
	elif os.path.isfile(pheno) != True:
	        sys.stderr.write('\nERROR: '+pheno+' not found\n')
	        sys.exit(1)

	#arguments
	Kpop = h5py.File(Kpop, 'r')
	pheno = h5py.File(pheno, 'r')
	hidden_k = int(hidden_k)
	snr_threshold = int(snr_threshold)
	fileout = h5py.File(fileout, 'w')

	#read Kpop file
	phenotype = pheno['phenotype/Ytransformed'][:] #catch warnings here
	kinship = Kpop['Kpop'][:] #catch warnings here	
	samples_rows = Kpop['Kpop/row_header/sample_ID'][:] #catch warnings here
	samples_columns = Kpop['Kpop/col_header/sample_ID'][:] #catch warnings here 
	#standardise kinship
	kinship /= sp.diag(kinship).mean()
	kinship += 1e-4 * sp.eye(kinship.shape[0])
	
	#select most variable genes
	selected_phenotypes = select_genes_on_SNR(phenotype, snr_threshold)
	print "\n{0} genes selected, featuring a signal-to-noise ratio >= {1} percentile\n".format(selected_phenotypes.shape[1],snr_threshold)

	#apply panama
	#cum_var = PC_varExplained(selected_phenotype)
	p = panama.PANAMA(Y=selected_phenotypes[:], Kpop = kinship[:])
	p.train(rank=hidden_k)
	Kpanama = p.get_Kpanama()
	Ktot = p.get_Ktot()
	vComp = p.get_varianceComps()
	print "\nVariance composition:\n{0}".format(vComp)
	print "\n Creating output file...\n"
	#writing into fileout	
	Ktot_shape=Ktot[:].shape
	dset=fileout.create_dataset('Ktot', Ktot_shape, dtype='float64')
	dset[...]=Ktot[:]
	dset1=fileout.create_dataset('Ktot/row_header/sample_ID', samples_rows.shape, dtype='S1000')
	dset1[...]=samples_rows
	dset2=fileout.create_dataset('Ktot/col_header/sample_ID', samples_columns.shape,dtype='S1000')
	dset2[...]=samples_columns

	fileout.close()
	print "\nCreating output file... done.\n"

	sys.exit(0)
