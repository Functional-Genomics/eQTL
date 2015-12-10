#!/usr/bin/env python

import sys
import scipy as SP
import h5py
import pdb
import copy
import warnings


class data():
	def __init__(self,geno,kinship,pheno,cov,cor,cormethod,window):
		"""load data file"""
		self.g	= h5py.File(geno,'r') #import geno data
		self.k = h5py.File(kinship,'r') # import kinship
		self.p = h5py.File(pheno,'r') #import pheno data
		self.c = h5py.File(cov,'r') #import covariates matrix
		self.correction = h5py.File(cor,'r') # import residuals from peer | Ktot from panama | Kpop iif no correction selected
		self.geneID = self.p['phenotype']['col_header']['phenotype_ID'][:] # ENSEMBL genes
		self.corrmeth = cormethod #string [ peer | panama | none ]
		self.window = float(window) #float value with nt window
	def getGeneIDs(self):
		""" get geneIDs """
		_chr = self.p['phenotype/chrom'][:] #upload vector of chr for each gene
		Iin = (_chr!='X')*(_chr!='Y')*(_chr!='MT') #here we are excluding these chromosomes from the analysis!
		rv = self.geneID[Iin]
		return rv 

	def getK(self,Isample=None,normalize=False):
		"""
		get Ktot/Kpop for sample specified in Isample
		"""
		if self.corrmeth == 'peer':
			RV = self.k['Kpop'][:]
			if normalize:
				RV /=RV.diagonal().mean()
		elif self.corrmeth == 'panama':
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

	def getGeneExpression(self,geneID,standardize=False):
		"""
		Get gene expression levels
		"""
		idx = self.geneID==geneID
		if self.corrmeth == 'peer':
			Y = self.correction['phenotype'][:,idx]
		elif self.corrmeth  == 'panama':
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
		

	def getGermlineExpr(self,geneID,standardize=False,Is=None,debug=False):
		"""
		get genotypes, chrom, pos. TODO1: change the code for trans analysis 
		"""
		w = self.window
		genePos = self.getGenePos(geneID)
		pos = self.g['genotype/col_header/pos'][:]
		chrom = self.g['genotype/col_header/chrom'][:]
		if w == 0:
			#pairwise comparison of each gene with all the SNPs
			X = self.g['genotype/matrix'][:]
			info = {}
			for key in self.g['genotype/col_header'].keys():
				info[key] = self.g['genotype/col_header'][key][:]
			return X, info
		else:
			Icis  = (chrom==float(genePos[0])) # force comparison on the same chromosome of the gene
			Icis *= (pos>float(genePos[1])-w) # select downstream cis SNPs
			Icis *= (pos<float(genePos[2])+w) # select upstream cis SNPs
			assert Icis.sum()>0, 'no cis intersection found'
			X = self.g['genotype/matrix'][:,Icis]
			info = {}
			for key in self.g['genotype/col_header'].keys():
				info[key] = self.g['genotype/col_header'][key][:][Icis] 
			return X, info



