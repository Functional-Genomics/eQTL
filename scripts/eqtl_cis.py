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

eqtl_cis.py <chr1.hdf5> <pheno.filtered.hdf5> <peer> <peer.hdf5> <Kpop.hdf5> <covariates.hdf5> <peer_cov> <cis_window> <n_perm> <nfolds> <fold_j> <outfilename>

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
#set the widget for the progressbar
bar = progressbar.ProgressBar(maxval=n_perm, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
#initialise the bar
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
		#initialize an array with one shape  = n_perms to store for each permutation the minimum pvalue per gene
		RV['pv0_min'] = SP.zeros(n_perm)
		SP.random.seed(0) #set random seed for each gene
		for perm_i in xrange(int(n_perm)):
			#print 'computing permutation # {0}'.format(perm_i)
			idx = SP.random.permutation(Xc.shape[0]) #take indexes
			Xc_perm = Xc[idx,:] #shuffle the samples of the genome matrix
			lmm_perm =run_lmm(booleanK,peer_cov,Xc_perm,Y,cov,K) #run the lmm model on permuted genotype
			pv_perm = lmm_perm.getPv() #get permuted pvalues
			pv0 = pv_perm[0,:] 
			RV['pv0_min'][perm_i] = pv0.min() #take the minimum pvalue of the list and populate the array
			if perm_i == 0:
				RV['lambda_perm'] = getLambda(pv_perm) #calculate lambda on the first permutation
			bar.update(perm_i + 1) #update the bar
		bar.finish()
	else:	
		print 'number of permutations is set = 1; empirical pvalues will not be computed.'
		lmm_perm =run_lmm(booleanK,peer_cov,Xc,Y,cov,K)

	#store results
	RV['pv'] = pv #record nominal pvalues
	if n_perm > 1:
		#compute how many MINIMUM permuted pvalues for each permutation are less than each nominal pvalue and store the value
		RV['pv_perm'] = SP.array([(sum(RV['pv0_min'][:]<pv[0,nominal]) for nominal in xrange(pv[:].shape[1])],dtype=float)
		RV['pv_perm'] += 1 # compute the empirical pvalues
		RV['pv_perm'] /= float(n_perm)+1 #compute the empirical pvalues
		RV['pv_perm'] = RV['pv_perm'].reshape((1,len(RV['pv_perm']))) #reshape
	else:
		RV['pv_perm'] = lmm_perm.getPv() #record just the permuted_pvalues after 1 permutation

	
	RV['qv'] = FDR.qvalues(pv) #multiple test correction for nominal pvalues
	RV['lambda'] = getLambda(pv) #get lambda for nominal pvalues
	if n_perm > 1:
		RV['lambda_empirical'] = getLambda(RV['pv_perm'])
	else:	
		RV['lambda_perm'] = getLambda(RV['pv_perm'])
		RV['lambda_empirical'] = (SP.empty((1,))).astype(str)
		RV['lambda_empirical'][:] = 'NA' # create an array with NA

	RV['beta'] = lmm.getBetaSNP() #get beta on nominal pvalues. TODO: should I get beta on permuted?
			
 	# add gene info
	for key in geno_info.keys():
	    RV[key] = geno_info[key]

	# export
	gene_group = fout.create_group(gene)
	dumpDictHdf5(RV,gene_group)
	print "gene kept %s"%gene

fout.close()

sys.exit(0)


