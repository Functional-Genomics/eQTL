#!/usr/bin/env python 

import sys

from utils import smartDumpDictHdf5
from utils import dumpDictHdf5
from utils import getLambda

import progressbar
import limix.modules.qtl as QTL
import limix.stats.fdr as FDR

import readline # because of an error with importing rpy2.robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

import scipy as SP
import scipy.linalg as LA
import os
import cPickle
import pdb
import time
import h5py

stats = importr('stats')

def usage():
        print '''
This script runs the eqtl analysis on a chunk of the gene expression matrix.

Usage:

eqtl_cis.py <chr1.hdf5> <pheno.filtered.hdf5> <peer> <peer.hdf5> <Kpop.hdf5> <covariates.hdf5> <use_kinship> <peer_cov> <cis_window> <n_perm> <change_beta_sign> <nfolds> <fold_j> <outfilename> <info_perm_file> <multiple_test_correction> <gene_cov.hdf5> 

peer_cov values = n | y [default=n] 
use_kinship = n | y 
gene_cov.hdf5 is optional [default=n]
multiple_test_correction= fdr | bonferroni [default=fdr]
'''

def run_lmm(use_kinship,peer_cov,Xc,Y,cov,K):
	'''this functions computes the lmm model using X (genotype), Y(phenotype). Optional: cov(matrix of covariance), K(kinship)'''
	if use_kinship=='y': #if use kinship, so K=K
		if peer_cov =='n': #if no covariates were used with peer then account for cov in the model
			sys.stderr.write('\nrunning lmm = QTL.test_lmm(Xc,Y,covs=cov)...')
			lmm = QTL.test_lmm(Xc,Y,covs=cov,K=K)
		else:
			sys.stderr.write('\nrunning lmm = QTL.test_lmm(Xc,Y,K=K)...')
			lmm = QTL.test_lmm(Xc,Y,K=K) # otherwise exclude covariates from the model since already used in peer
	else: #if not use kinship
		if peer_cov =='n': #if no cov where used with peer then
			sys.stderr.write('\nrunning lmm= QTL.test_lmm(Xc,Y,covs=cov,K=SP.eye(Xc.shape[0]))...')
			lmm = QTL.test_lmm(Xc,Y,covs=cov, K=SP.eye(Xc.shape[0])) #use cov in the model
		else:
			sys.stderr.write('\nrunning lmm = QTL.test_lmm(Xc,Y,K=SP.eye(Xc.shape[0]))...')
			lmm = QTL.test_lmm(Xc,Y, K=SP.eye(Xc.shape[0])) #exclude cov in the model since already used by peer
	return lmm

if len(sys.argv[1:]) < 16:
	sys.stderr.write('ERROR: missing parameters\n')
	usage()
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
	sys.stderr.write('\nERROR: Please use either y or n for beta sign\n')
	sys.exit(1)
elif use_kinship != 'y' and use_kinship != 'n':
        usage()
        sys.stderr.write('\nERROR: Please use either y or n for use_kinship\n')
        sys.exit(1)

#open outfile 
fout  = sys.argv[14]  #%d_%.3d.hdf5'%(nfolds,fold_j) #this should take as argument a name like nfolds_j.hdf5
info_out = sys.argv[15] #hdf5 file to store permuted pvalues and inflation factors
multiple_test_correction = sys.argv[16] #multiple test correction method chosen

if len(sys.argv[1:]) > 17:
	usage()
	sys.stderr.write("ERROR: Too many arguments\n")
	sys.exit(1)


if len(sys.argv[1:]) == 17:
	gene_cov = h5py.File(sys.argv[17],'r')
	gene_cov_matrix = gene_cov['covariates'][:] #get matrix of covariates
	gene_cov_genes = gene_cov['col_header/phenotype_ID'][:] #get phenotype names from gene covariates
	g_cov_used = 'y'
else:
	g_cov_used = 'n'

sys.stderr.write('\nParameters set in the association analysis:\ngenotype={0}; phenotype={1}; correction_metho={2}; correction_method_file={3}; Kinship_file={4}; covariates_file={5}; use_kinship={6}; peer_covariates={7}; cis_window={8}; number_permutation={9}; change_beta_sign={10}; n_folds={11}; fold_j={12}; outfile={13}; gene_covariates_used={14}'.format(geno,pheno,cm,cm_hdf5,kinship,cov_hdf5,use_kinship,peer_cov,window,n_perm,change_beta_sign,nfolds,fold_j,fout,g_cov_used))

#open the outfile
fout = h5py.File(fout,'w')
info_out = h5py.File(info_out,'w')
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
if use_kinship=='y':
	#check if Kpop or Ktot contains Nan
	booleanK=SP.isnan(K)
	if True in booleanK:
		sys.stderr.write('\nWARNING: kinship matrix contains NA\n')

#set the widget for the progressbar
bar = progressbar.ProgressBar(maxval=n_perm, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
#execute analysis for each gene in chunk j
for gene in genes:

	print "\n.. gene %s"%gene

	#1. get geno and pheno data
	Y = data.getGeneExpression(gene,standardize=False)
	#try:
	Xc,geno_info = data.gene_SNP_pair(gene)
	if len(Xc)==0:
		print "-"
		continue
#	except:
#	    e = sys.exc_info()[0]
#	    print "-"
#	    continue
	
	print "+"
	
	if g_cov_used == 'y': #if gene covariates are used
		gene_split = gene.split('.')[0]
		#i = SP.in1d(gene_cov_genes[:],gene)
		gene_cov_genes_split = SP.array(map(lambda x:x.split('.')[0],gene_cov_genes.tolist()))
		i = SP.in1d(gene_cov_genes_split[:],gene_split)
		igene_cov_matrix = gene_cov_matrix[:,i]
		g_cov = SP.hstack((cov[:],igene_cov_matrix[:]))
	else:
		g_cov = cov[:]
	#run the linear mixed model to get nominal pvalues
	lmm = run_lmm(use_kinship,peer_cov,Xc,Y,g_cov,K)
	pv = lmm.getPv() #store nominal pv
	RV={} #open empty dictionary to store results
	INFO={} #open empty dictionary to store permuted pvalues

	#permutation	
	if n_perm > 100:
		print '\nnumber of permutations is set > 100; empirical pvalues will be computed.'
		#initialize the bar
		bar.start()
		#initialize an array with one shape  = n_perms to store for each permutation the minimum pvalue per gene
		r = 0 #initialize fraction of permuted pvalues <= minimum nominal pval
		SP.random.seed(0) #set random seed for each gene
		for perm_i in xrange(int(n_perm)):
			#print 'computing permutation # {0}'.format(perm_i)
			idx = SP.random.permutation(Xc.shape[0]) #take indexes
			Xc_perm = Xc[:][idx,:] #shuffle the samples of the genome matrix
			lmm_perm =run_lmm(use_kinship,peer_cov,Xc_perm,Y,g_cov,K) #run the lmm model on permuted genotype
			pv_perm = lmm_perm.getPv() #get permuted pvalues
			pv0_min= pv_perm[0,:].min() #take the minimum pvalue of the list 
			if pv0_min <= pv[:].min(): #if minimum permuted pval is less or equal than the nominal increase r
				r+=1
			if perm_i == 0:
				RV['lambda_perm'] = getLambda(pv_perm) #calculate lambda on the first permutation
			bar.update(perm_i + 1) #update the bar
		bar.finish()
	else:	
		print '\nnumber of permutations is set <=100; empirical pvalues will not be computed.'
		SP.random.seed(0)
		list_perm_pvalues = []
		list_perm_lambdas = []
		for perm_i in xrange(int(n_perm)):
			idx = SP.random.permutation(Xc.shape[0])
			Xc_perm = Xc[:][idx,:]
			lmm_perm =run_lmm(use_kinship,peer_cov,Xc_perm,Y,g_cov,K)
			perm_pv = lmm_perm.getPv()
			if perm_i == 0:
				perm_pv_tmp_array = perm_pv
			else:
				perm_pv_tmp_array = SP.concatenate((perm_pv_tmp_array,perm_pv),axis=0)
			pv0_min = perm_pv[0,:].min()
			list_perm_pvalues.append(pv0_min)
			list_perm_lambdas.append(getLambda(perm_pv[:])[0][0])
		list_perm_pvalues = SP.array(list_perm_pvalues)
		min_perm_pvalue = list_perm_pvalues.min()
		list_perm_lambdas = SP.array(list_perm_lambdas)
		mean_perm_lambda = list_perm_lambdas.mean()
		std_perm_lambda = list_perm_lambdas.std()
		INFO['all_pv_perm'] = perm_pv_tmp_array[:]
		INFO['all_lambda_perm'] = list_perm_lambdas[:]
		
	#store results
	RV['pv'] = pv #record nominal pvalues
	if n_perm > 100:
		#compute how many MINIMUM permuted pvalues for each permutation are less than the minimum nominal pvalue and store the value
		RV['pv_perm'] = SP.array([((r+1)/(float(n_perm)+1))],dtype=float) #one value per gene
		RV['pv_perm'] = RV['pv_perm'].reshape((1,len(RV['pv_perm']))) #reshape an array of (1,1)
	else: #if n_perm <= 100 
		#perm_pv= lmm_perm.getPv() #record the minimum permuted_pvalues after 1 permutation
		#RV['pv_perm'] = SP.array([perm_pv[:].min()])
		RV['pv_perm'] = SP.array([min_perm_pvalue])

	if multiple_test_correction == 'fdr': #default
		RV['qv'] = FDR.qvalues(pv) #multiple test correction for nominal pvalues with B-H
	else:
		RV['qv'] = SP.empty(pv.shape,dtype=float)
		tmp_qv = SP.array(stats.p_adjust(FloatVector(pv[0].tolist()),method = 'bonferroni',n=float(pv.shape[1]))) #multiple test correction for nominal pvalues with Bonferroni 
		RV['qv'][0,:]=tmp_qv
	RV['lambda'] = getLambda(pv) #get lambda for nominal pvalues
	if n_perm <= 100:
		#RV['lambda_perm'] = getLambda(perm_pv[:]) #calculate lambda on the the permutation
		RV['lambda_perm'] = SP.array([mean_perm_lambda,std_perm_lambda])

	if change_beta_sign == 'y':
		RV['beta'] = -lmm.getBetaSNP() #get beta on nominal pvalues and change sign: this gives always the effect size computed between the ALT with respect to the REF allele as defined in the VCF
	else:
		RV['beta'] = lmm.getBetaSNP() #get beta on nominal pvalues. This gives always the effect size computed between the REF with respect to the ALT allele as defined in the VCF
			
 	# add gene info
	for key in geno_info.keys():
	    RV[key] = geno_info[key]

	# export
	gene_group = fout.create_group(gene)
	gene_group2 = info_out.create_group(gene)
	dumpDictHdf5(RV,gene_group)
	dumpDictHdf5(INFO,gene_group2)
	print "\ngene kept %s"%gene

fout.close()

sys.exit(0)


