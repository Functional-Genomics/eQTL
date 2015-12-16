#!/usr/bin/env python 

import sys

from utils import smartDumpDictHdf5
from utils import dumpDictHdf5
from utils import getLambda

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

if len(sys.argv[1:]) < 12:
	usage()
	sys.stderr.write('ERROR: missing parameters\n')
	sys.exit(1)

#read args
geno,pheno,cm,cm_hdf5,kinship,cov_hdf5,peer_cov,window,n_perm = sys.argv[1:10]
window=float(window)
#populate dictionary with all the data needed for eqtl analysis
#rom eqtlsettings import read_args as ra
#CFG,correction_method = ra(geno = sys.argv[1], pheno=sys.argv[2], correction_method = sys.argv[3], hdf5_correction =sys.argv[4], Kpop = sys.argv[5], covariates = sys.argv[6])


#take nfold and j to name the out file for each j
nfolds = int(sys.argv[10])
fold_j = int(sys.argv[11])

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

#execute analysis for each gene in chunk j
for gene in genes:

	print ".. gene %s"%gene

	#1. get geno and pheno data
	Y = data.getGeneExpression(gene,standardize=False)
	try:
	    Xc,geno_info = data.getGermlineExpr(gene)
	except:
	    e = sys.exc_info()[0]
	    print "...excluding gene %s %s" %(gene,e) 
	    continue
	#permutation

	for perm_i in SP.arange(int(n_perm)):
		SP.random.seed(perm_i)
		idx = SP.random.permutation(Xc.shape[0])
		Xc_perm = Xc[idx,:]

	#check if Kpop or Ktot contains Nan
	booleanK=SP.isnan(K)
	if True in booleanK:
		if peer_cov =='n': #if no covariates were used with peer then account for cov in the model
			lmm = QTL.test_lmm(Xc,Y,covs=cov)
			lmm_perm = QTL.test_lmm(Xc_perm,Y,covs=cov)
		else:
			lmm = QTL.test_lmm(Xc,Y) # otherwise exclude covariates from the model since already used in peer
			lmm_perm = QTL.test_lmm(Xc_perm,Y)
	else:
		if peer_cov =='n': #if no cov where used with peer then
			lmm = QTL.test_lmm(Xc,Y,covs=cov,K=K) #use cov in the model
			lmm_perm = QTL.test_lmm(Xc_perm,Y,covs=cov,K=K)
		else:
			lmm = QTL.test_lmm(Xc,Y,K=K) #exclude cov in the model since already used by peer
			lmm_perm = QTL.test_lmm(Xc_perm,Y,K=K)

	# run the linear mixed model
	pv=lmm.getPv()
	pv_perm = lmm_perm.getPv()
	RV = {}
	RV['pv'] = pv
	RV['pv_perm'] = pv_perm
	RV['qv'] = FDR.qvalues(pv)
	RV['lambda'] = getLambda(pv)
	RV['lambda_perm'] = getLambda(pv_perm)
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


