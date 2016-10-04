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

eqtl_trans.py <allchr.hdf5> <pheno.filtered.hdf5> <peer> <peer.hdf5> <Kpop.hdf5> <covariates.hdf5> <use_kinship> <peer_cov> <cis_window> <n_perm> <change_beta_sign> <nfolds> <fold_j> <outfilename>

peer_cov values = n | y [default=n] 
use_kinship = n| y
'''

def run_lmm(use_kinship,peer_cov,Xc,Y,cov,K):
        '''this functions computes the lmm model using X (genotype), Y(phenotype). Optional: cov(matrix of covariance), K(kinship)'''
	if use_kinship=='y':
		if peer_cov =='n': #if no covariates were used with peer then account for cov in the model
			sys.stderr.write('\nrunning lmm = QTL.test_lmm(Xc,Y,covs=cov,K=K)...')
			lmm = QTL.test_lmm(Xc,Y,covs=cov,K=K)
		else:
			sys.stderr.write('\nrunning lmm = QTL.test_lmm(Xc,Y,K=K)...')
			lmm = QTL.test_lmm(Xc,Y,K=K) # otherwise exclude covariates from the model since already used in peer
	else: #no kinship
		if peer_cov =='n': #if no cov where used with peer then
			sys.stderr.write('\nrunning lmm= QTL.test_lmm(Xc,Y,covs=cov,K=SP.eye(Xc.shape[0]))...')
			lmm = QTL.test_lmm(Xc,Y,covs=cov,K=SP.eye(Xc.shape[0])) #include covariates in the model
		else:
			sys.stderr.write('\nrunning lmm = QTL.test_lmm(Xc,Y,K=SP.eye(Xc.shape[0]))...')
			lmm = QTL.test_lmm(Xc,Y,K=SP.eye(Xc.shape[0])) #exclude cov in the model since already used by peer
	return lmm


if len(sys.argv[1:]) < 14:
	usage()
	sys.stderr.write('ERROR: missing parameters\n')
	sys.exit(1)

#read args
geno,pheno,cm,cm_hdf5,kinship,cov_hdf5,use_kinship,peer_cov,window,n_perm,change_beta_sign = sys.argv[1:12]
window=float(window)
n_perm=int(n_perm)

#populate dictionary with all the data needed for eqtl analysis

#take nfold and j to name the out file for each j
nfolds = int(sys.argv[12])
fold_j = int(sys.argv[13])


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
	sys.stderr.write('\nERROR: Please use either y or n for beta sign option\n')
	sys.exit(1)
elif use_kinship !='y' and use_kinship !='n':
	usage()
	sys.stderr.write('\nERROR: Please use either y or n for use_kinship\n')
	sys.exit(1)

#open outfile 
fout  = h5py.File(sys.argv[14],'w')  #%d_%.3d.hdf5'%(nfolds,fold_j) #this should take as argument a name like nfolds_j.hdf5

sys.stderr.write('\nParameters set in the association analysis:\ngenotype={0}; phenotype={1}; correction_metho={2}; correction_method_file={3}; Kinship_file={4}; covariates_file={5}; use_kinship={6}; peer_covariates={7}; cis_window={8}; number_permutation={9}; change_beta_sign={10}; n_folds={11}; fold_j={12}; outfile={13}'.format(geno,pheno,cm,cm_hdf5,kinship,cov_hdf5,use_kinship,peer_cov,window,n_perm,change_beta_sign,nfolds,fold_j,fout))

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
#output warning if NAN in Kinship
if use_kinship=='y':
	booleanK=SP.isnan(K)
	if True in booleanK:
		sys.stderr.write('\nWARNING: kinship matrix contains NA\n')

#set progress bar
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
#		print "...excluding gene %s %s" %(gene,e)

        print "+"

 	#run the linear mixed model to get nominal pvalues
	lmm = run_lmm(use_kinship,peer_cov,Xc,Y,cov,K)
	pv = lmm.getPv() #store nominal pv
	RV={} #open empty dictionary to store results
	
	#permutation
	r = 0
	if n_perm > 1:
		print 'number of permutations is set > 1; empirical pvalues will be computed.'
		#initilize the bar
		bar.start()
		SP.random.seed(0) #set random seed for each gene
		for perm_i in xrange(int(n_perm)):
			idx = SP.random.permutation(Xc.shape[0]) #take indexes
			Xc_perm = Xc[idx,:] #shuffle the samples of the genome matrix
			lmm_perm =run_lmm(use_kinship,peer_cov,Xc_perm,Y,cov,K) #run the lmm model on permuted genotype
			pv_perm = lmm_perm.getPv()# get pvalues 
			pv0_min = pv_perm[0,:].min() #get minimum permuted pvalue
			if pv0_min <= pv[:].min():
				r+=1 #increment r
			if perm_i == 0: #get lambda at the first permutation
				RV['lambda_perm'] = getLambda(pv_perm)
			bar.update(perm_i+1)
		bar.finish()
	else:
		print 'number of permutations is set = 1; empirical pvalues will not be computed.'
		idx = SP.random.permutation(Xc.shape[0]) #take indexes
		Xc_perm = Xc[idx,:] #shuffle the samples of the genome matrix
		lmm_perm =run_lmm(use_kinship,peer_cov,Xc_perm,Y,cov,K)

	# run the linear mixed model
	RV['pv'] = pv
	#compute  qvalues with Benjamini - Hochberg
	RV['qv'] = (SP.empty((1,pv.shape[1]))).astype(str)
	RV['qv'][:] = 'NA' #fill the local adjusted pvalues with NAN since multiple test correction will be applied globally
	#compute empirical pvalues if n_perm > 1
	if n_perm > 1:
		#compute how many MINIMUM permuted pvalues for each permutation are greater than each nominal pvalue and store the value
		RV['pv_perm'] = SP.array([((r+1)/(float(n_perm)+1))],dtype=float) #generate array with one empirical pvaluee per gene
		RV['pv_perm'] = RV['pv_perm'].reshape((1,len(RV['pv_perm']))) #reshape
	else:
		perm_pv = lmm_perm.getPv() #do not get empirical pvalue but just permuted ones
		perm_pv_min = perm_pv[:].min()
		RV['pv_perm'] = SP.array([perm_pv_min])

	
	#compute lambda
	RV['lambda'] = getLambda(pv) # compute lambda on nominal pvalues
	if n_perm ==  1: #no empirical pvalues
		RV['lambda_perm'] = getLambda(perm_pv[:]) #get lambda for permuted pvalues

	#compute beta
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


