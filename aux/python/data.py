#!/usr/bin/env python

import sys
import scipy as SP
import h5py
import pdb
import copy
import warnings

#from eqtlsettings import read_args

#CFG,correction_method=read_args(data_argv[0], data_argv[1], data_argv[2], data_argv[3], data_argv[4], data_argv[5])


class data():
	def __init__(self):
		self.load()

	def load(self):
		"""load data file:
		cache_genotype: load genotypes fully into memory (False)
		"""
		self.g	= h5py.File(CFG['data']['geno'],'r') #import geno data
		self.k = h5py.File(CFG['data']['kinship'],'r') # import kinship
		self.p = h5py.File(CFG['data']['pheno'],'r') #import pheno data
		self.c = h5py.File(CFG['data']['covariates'],'r') #import covariates matrix
		self.correction = h5py.File(CFG['data']['correction'],'r') # import residuals from peer | Ktot from panama | Kpop iif no correction selected
		self.geneID = self.p['phenotype']['col_header']['phenotype_ID'][:] # ENSEMBL genes

	def getGeneIDs(self):
		""" get geneIDs """
		_chr = self.p['phenotype/chrom'][:] #upload vector of chr for each gene
		Iin = (_chr!='X')*(_chr!='Y')*(_chr!='MT') #here we are excluding these chromosomes from the analysis!
		rv = self.geneID[Iin]
		return rv 

	def getK(self,Isample=None,correction_method = correction_method,normalize=False):
		"""
		get Ktot/Kpop for sample specified in Isample
		"""
		if correction_method == 'peer':
			RV = self.k['Kpop'][:]
			if normalize:
				RV /=RV.diagonal().mean()
		elif correction_method == 'panama':
			RV = self.correction['Ktot'][:]
		else :
			RV =self.k['Kpop'][:]
			if normalize:
				RV /=RV.diagonal().mean()
		if Isample!=None:
			RV = RV[Isample,:][:,Isample]
		return RV

	def getCovariates(self):
		"""
		get covariates
		"""
		rv = self.c['covariates'][:]
		#col = 1.*(rv.sum(1)==0)[:,SP.newaxis]
		#rv = SP.concatenate([rv,col],1)
		return rv

	def getGeneExpression(self,geneID,correction_method = correction_method,standardize=False):
		"""
		Get gene expression levels
		"""
		idx = self.geneID==geneID
		if correction_method == 'peer':
			Y = self.correction['phenotype'][:,idx]
		elif correction_method == 'panama':
			Y = self.p['phenotype/Ytransformed'][:,idx]
		else:
			Y = self.p['phenotype/Ytransformed'][:,idx] #when none is selected
		if standardize: #at the moment this standardization is not optional in the pipeline. TODO: make it optional
			Y-=Y.mean(0)
			Y/=Y.std(0)
		return Y

	def getGenePos(self,geneID):
		"""
		get position of the gene
		"""
		idx = SP.where(self.geneID==geneID)[0][0]
		gene_chrom = self.p['phenotype/chrom'][idx] 
		gene_start = float(self.p['phenotype/start'][idx]) 
		gene_end = float(self.p['phenotype/end'][idx])
		rv = SP.array([gene_chrom,gene_start,gene_end])
		return rv
		

	def getGermlineExpr(self,geneID,cis_window=1000000.0,standardize=False,Is=None,debug=False):
		"""
		get genotypes, chrom, pos. TODO: cis window is not optional at the moment; TODO1: trans analysis is not optional! 
		"""
		genePos = self.getGenePos(geneID)
		pos = self.g['genotype/col_header/pos'][:]
		chrom = self.g['genotype/col_header/chrom'][:]
		Icis  = (chrom==float(genePos[0]))
		Icis *= (pos>float(genePos[1])-cis_window)
		Icis *= (pos<float(genePos[2])+cis_window)
		assert Icis.sum()>0, 'no cis intersection found'
		X = self.g['genotype/matrix'][:,Icis]
		info = {}
		for key in self.g['genotype/col_header'].keys():
			info[key] = self.g['genotype/col_header'][key][:][Icis] 
		return X, info



