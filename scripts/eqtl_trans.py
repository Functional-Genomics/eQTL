#!/usr/bin/env python 

import sys

from utils import smartDumpDictHdf5
from utils import dumpDictHdf5
from utils import getLambda

import limix.modules.qtl as QTL
import limix.stats.fdr as FDR
import progressbar
import scipy as SP
import scipy.linalg as LA
import os
import cPickle
import pdb
import time
import h5py

def usage():
        print '''
This script runs the eqtl analysis on a chunk of the gene expression matrix VS all the SNPs.

Usage:

eqtl_trans.py <allchr.hdf5> <pheno.filtered.hdf5> <peer> <peer.hdf5> <Kpop.hdf5> <covariates.hdf5> <peer_cov> <cis_window> <n_perm> <nfolds> <fold_j> <outfilename>

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


if len(sys.argv[1:]) < 12:
	usage()
	sys.stderr.write('ERROR: missing parameters\n')
	sys.exit(1)

#read args
geno,pheno,cm,cm_hdf5,kinship,cov_hdf5,peer_cov,window,n_perm = sys.argv[1:10]
window=float(window)
n_perm=int(n_perm)

#populate dictionary with all the data needed for eqtl analysis
#rom eqtlsettings import read_args as ra
#CFG,correction_method = ra(geno = sys.argv[1], pheno=sys.argv[2], correction_method = sys.argv[3], hdf5_correction =sys.argv[4], Kpop = sys.argv[5], covariates = sys.argv[6])

#take nfold and j to name the out file for each j
nfolds = int(sys.argv[10])
fold_j = int(sys.argv[11])


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

#open outfile 
fout  = h5py.File(sys.argv[12],'w')  #%d_%.3d.hdf5'%(nfolds,fold_j) #this should take as argument a name like nfolds_j.hdf5
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
#set progress bar
bar = progressbar.ProgressBar(maxval=n_perm, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
bar.start()
#execute analysis for each gene in chunk j
for gene in genes:

	print ".. gene %s"%gene

	#1. get geno and pheno data
	Y = data.getGeneExpression(gene,standardize=False)
	try:
		Xc,geno_info = data.gene_SNP_pair(gene)
	except:
		e = sys.exc_info()[0]
		print "...excluding gene %s %s" %(gene,e) 
		continue

        #check if Kpop or Ktot contains Nan
	booleanK=SP.isnan(K)
 	#run the linear mixed model to get nominal pvalues
	lmm = run_lmm(booleanK,peer_cov,Xc,Y,cov,K)
	pv = lmm.getPv() #store nominal pv
	RV={} #open empty dictionary to store results
	
	#permutation
	if n_perm > 1:
		print 'number of permutations is set > 1; empirical pvalues will be computed.'
		#initialize an array of 0 with shape of n_perms X number of SNPs in the cis window
		#RV['pv0'] = SP.zeros((n_perm,pv.shape[1]))
		#initialize an array with one shape  = n_perms to store for each permutation the minimum pvalue per gene
		RV['pv0_min'] = SP.zeros(n_perm)
		SP.random.seed(0) #set random seed for each gene
		for perm_i in xrange(int(n_perm)):
			idx = SP.random.permutation(Xc.shape[0]) #take indexes
			Xc_perm = Xc[idx,:] #shuffle the samples of the genome matrix
			lmm_perm =run_lmm(booleanK,peer_cov,Xc,Y,cov,K) #run the lmm model on permuted genotype
			pv_perm = lmm_perm.getPv()# get pvalues 
			pv0 = pv_perm[0,:] #get permuted pvalue
			#RV['pv0'][perm_i,:] = pv0 #populate the array with list of permuted pvalues in each row
			RV['pv0_min'][perm_i] = pv0.min() #take the minimum pvalue of the list and populate the array
			if perm_i == 0: #get lambda at the first permutation
				RV['lambda_perm'] = getLambda(pv_perm)
			bar.update(perm_i+1)
		bar.finish()
	else:
		print 'number of permutations is set = 1; empirical pvalues will not be computed.'
		lmm_perm =run_lmm(booleanK,peer_cov,Xc,Y,cov,K)

	# run the linear mixed model
	RV['pv'] = pv
	#compute  qvalues with Benjamini - Hochberg
	RV['qv'] = FDR.qvalues(pv) 
	#compute empirical pvalues if n_perm > 1
	if n_perm > 1:
		#compute how many MINIMUM permuted pvalues for each permutation are less than each nominal pvalue and store the value
		RV['pv_perm'] = SP.array([(RV['pv0_min']<pv[0,nominal]).sum() for nominal in xrange(pv.shape[1])],dtype=float)
		RV['pv_perm'] += 1 # compute the empirical pvalues
		RV['pv_perm'] /= float(n_perm) #compute the empirical pvalues
		RV['pv_perm'] = RV['pv_perm'].reshape((1,len(RV['pv_perm']))) #reshape
	else:
		RV['pv_perm'] = lmm_perm.getPv() #do not get empirical pvalue but just permuted ones

	
	#compute lambda
	RV['lambda'] = getLambda(pv) # compute lambda on nominal pvalues
	if n_perm ==  1: #no empirical pvalues
		RV['lambda_perm'] = getLambda(RV['pv_perm']) #get lambda for permuted pvalues
		RV['lambda_empirical'] = SP.empty((1,))
		RV['lambda_empirical'][:] = SP.NAN # create an array with NaN 
	else:	
		RV['lambda_empirical'] = getLambda(RV['pv_perm'])

	#compute beta
	RV['beta'] = lmm.getBetaSNP()
	
	# add gene info
	for key in geno_info.keys():
	    RV[key] = geno_info[key]

	# export
	gene_group = fout.create_group(gene)
	dumpDictHdf5(RV,gene_group)
	print "gene kept %s"%gene

fout.close()

sys.exit(0)


