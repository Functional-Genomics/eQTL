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

eqtl_cis.py <chr1.hdf5> <pheno.filtered.hdf5> <peer> <peer.hdf5> <Kpop.hdf5> <covariates.hdf5> <nfolds> <fold_j> <outfilename> '''

if len(sys.argv[1:]) < 9:
	sys.stderr.write('ERROR: missing parameters\n')
	usage()
	sys.exit(1)


#populate dictionary with all the data needed for eqtl analysis
from eqtlsettings import read_args as ra
CFG,correction_method = ra(geno = sys.argv[1], pheno=sys.argv[2], correction_method = sys.argv[3], hdf_correction =sys.argv[4], Kpop = sys.argv[5], covariates = sys.argv[6])

#take nfold and j to name the out file for each j
nfolds = int(sys.argv[7])
fold_j = int(sys.argv[8])

#load data
import data as DATA
#open outfile 
fout  = h5py.File(sys.argv[9],'w')  #%d_%.3d.hdf5'%(nfolds,fold_j) #this should take as argument a name like nfolds_j.hdf5
# load data 
data  = DATA.data()
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
	Y = data.getGeneExpression(gene,correction_method = correction_method, standardize=False)
	try:
	    Xc,geno_info = data.getGermlineExpr(gene,cis_window=1000000.0)
	except:
	    e = sys.exc_info()[0]
	    print "...excluding gene %s %s" %(gene,e) 
	    continue

	# run the linear mixed model
	lmm = QTL.test_lmm(Xc,Y,covs=cov,K=K)
	pv=lmm.getPv()
	RV = {}
	RV['pv'] = pv
	RV['qv'] = FDR.qvalues(pv)
	RV['lambda'] = getLambda(pv)
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


