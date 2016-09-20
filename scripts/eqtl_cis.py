#!/usr/bin/env python 

import sys

from utils import smartDumpDictHdf5
from utils import dumpDictHdf5
from utils import getLambda

import progressbar
import limix.modules.qtl as QTL
import limix.stats.fdr as FDR

import scipy as SP
import scipy.linalg as LA
import os
import cPickle
import pdb
import time
import h5py

def usage():
        print '''
This script runs the eqtl analysis on a chunk of the gene expression matrix.

Usage:

eqtl_cis.py <chr1.hdf5> <pheno.filtered.hdf5> <peer> <peer.hdf5> <Kpop.hdf5> <covariates.hdf5> <peer_cov> <cis_window> <n_perm> <change_beta_sign> <nfolds> <fold_j> <outfilename>

peer_cov values = n | y [default=n] '''

def run_lmm(booleanK,peer_cov,Xc,Y,cov,K):
	'''this functions computes the lmm model using X (genotype), Y(phenotype). Optional: cov(matrix of covariance), K(kinship)'''
	if True in booleanK:
		if peer_cov =='n': #if no covariates were used with peer then account for cov in the model
			lmm = QTL.test_lmm(Xc,Y,covs=cov)
		else:
			lmm = QTL.test_lmm(Xc,Y) # otherwise exclude covariates from the model since already used in peer
	else:
		if peer_cov =='n': #if no cov where used with peer then
			lmm = QTL.test_lmm(Xc,Y,covs=cov,K=K) #use cov in the model
		else:
			lmm = QTL.test_lmm(Xc,Y,K=K) #exclude cov in the model since already used by peer
	return lmm

if len(sys.argv[1:]) < 13:
	usage()
	sys.stderr.write('ERROR: missing parameters\n')
	sys.exit(1)

#read args
geno,pheno,cm,cm_hdf5,kinship,cov_hdf5,peer_cov,window,n_perm,change_beta_sign = sys.argv[1:11]
window=float(window)
n_perm=int(n_perm)
#populate dictionary with all the data needed for eqtl analysis

#take nfold and j to name the out file for each j

nfolds = int(sys.argv[11])
fold_j = int(sys.argv[12])

if type(nfolds) != (int):
	usage()
	sys.stderr.write('\nERROR: Please use integer for nfolds\n')
	sys.exit(1)
elif type(fold_j) != (int):
	usage()
	sys.stderr.write('\nERROR: Please use integer for fold_j\n')
	sys.exit(1)
elif type(window) != (float) and type(window) != (int):
	usage()
	sys.stderr.write('\nERROR: Please use integer for window\n')
	sys.exit(1)
elif change_beta_sign != 'y' and change_beta_sign != 'n':
	usage()
	sys.stderr.write('\nERROR: Please use either y or n for beta sign\n')
	sys.exit(1)

#open outfile 
fout  = h5py.File(sys.argv[13],'w')  #%d_%.3d.hdf5'%(nfolds,fold_j) #this should take as argument a name like nfolds_j.hdf5
# load data 
import data as DATA
data  = DATA.data(geno,kinship,pheno,cov_hdf5,cm_hdf5,cm,window)
#get kinship
K  = data.getK(normalize=False) #at the moment normalisation is not optional. Kpop/Ktot will be always normalised
#get number of samples
N = K.shape[0]
#get covariates
cov   = data.getCovariates()
#get genes
genes = data.getGeneIDs()
#get number of genes
n_genes = genes.shape[0]
#split the gene expression matrix in chunks equal to nfolds.
Icv = SP.floor(nfolds*SP.arange(n_genes)/n_genes)
#take the j chunk
I = Icv==fold_j
#grab the genes in chunk j
genes = list(genes[I])
#check if Kpop or Ktot contains Nan
booleanK=SP.isnan(K)
#set the widget for the progressbar
bar = progressbar.ProgressBar(maxval=n_perm, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
#execute analysis for each gene in chunk j
for gene in genes:

	print ".. gene %s"%gene

	#1. get geno and pheno data
	Y = data.getGeneExpression(gene,standardize=False)
	try:
	    Xc,geno_info = data.gene_SNP_pair(gene)
	except:
	    e = sys.exc_info()[0]
	    print "-"
	    continue
	
	print "+"
	
	
	#run the linear mixed model to get nominal pvalues
	lmm = run_lmm(booleanK,peer_cov,Xc,Y,cov,K)
	pv = lmm.getPv() #store nominal pv
	RV={} #open empty dictionary to store results

	#permutation	
	if n_perm > 1:
		print 'number of permutations is set > 1; empirical pvalues will be computed.'
		#initialize the bar
		bar.start()
		#initialize an array with one shape  = n_perms to store for each permutation the minimum pvalue per gene
		r = 0 #initialize fraction of permuted pvalues <= minimum nominal pval
		SP.random.seed(0) #set random seed for each gene
		for perm_i in xrange(int(n_perm)):
			#print 'computing permutation # {0}'.format(perm_i)
			idx = SP.random.permutation(Xc.shape[0]) #take indexes
			Xc_perm = Xc[:][idx,:] #shuffle the samples of the genome matrix
			lmm_perm =run_lmm(booleanK,peer_cov,Xc_perm,Y,cov,K) #run the lmm model on permuted genotype
			pv_perm = lmm_perm.getPv() #get permuted pvalues
			pv0_min= pv_perm[0,:].min() #take the minimum pvalue of the list 
			if pv0_min <= pv[:].min(): #if minimum permuted pval is less or equal than the nominal increase r
				r+=1
			if perm_i == 0:
				RV['lambda_perm'] = getLambda(pv_perm) #calculate lambda on the first permutation
			bar.update(perm_i + 1) #update the bar
		bar.finish()
	else:	
		print 'number of permutations is set = 1; empirical pvalues will not be computed.'
		idx = SP.random.permutation(Xc.shape[0])
		Xc_perm = Xc[:][idx,:]
		lmm_perm =run_lmm(booleanK,peer_cov,Xc_perm,Y,cov,K)

	#store results
	RV['pv'] = pv #record nominal pvalues
	if n_perm > 1:
		#compute how many MINIMUM permuted pvalues for each permutation are less than the minimum nominal pvalue and store the value
		RV['pv_perm'] = SP.array([((r+1)/(float(n_perm)+1))],dtype=float) #one value per gene
		RV['pv_perm'] = RV['pv_perm'].reshape((1,len(RV['pv_perm']))) #reshape an array of (1,1)
	else: #if n_perm ==1 
		perm_pv= lmm_perm.getPv() #record the minimum permuted_pvalues after 1 permutation
		RV['pv_perm'] = SP.array([perm_pv[:].min()])

	
	RV['qv'] = FDR.qvalues(pv) #multiple test correction for nominal pvalues
	RV['lambda'] = getLambda(pv) #get lambda for nominal pvalues
	if n_perm == 1:
		RV['lambda_perm'] = getLambda(perm_pv[:]) #calculate lambda on the the permutation

	if change_beta_sign == 'y':
		RV['beta'] = -lmm.getBetaSNP() #get beta on nominal pvalues and change sign: this gives always the effect size computed between the ALT with respect to the REF allele as defined in the VCF
	else:
		RV['beta'] = lmm.getBetaSNP() #get beta on nominal pvalues. This gives always the effect size computed between the REF with respect to the ALT allele as defined in the VCF
			
 	# add gene info
	for key in geno_info.keys():
	    RV[key] = geno_info[key]

	# export
	gene_group = fout.create_group(gene)
	dumpDictHdf5(RV,gene_group)
	print "gene kept %s"%gene

fout.close()

sys.exit(0)


