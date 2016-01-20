#!/usr/bin/env python

import sys
import limix.stats.fdr as FDR
from utils import smartAppend
from utils import smartDumpDictHdf5
import scipy as SP
import h5py
import os
import pdb
import glob
import cPickle


def usage():
	print '''
This script generates a final summary file from all the hdf5 files generated by the eqtl_cis/trans.py script. 
Usage:

eqtl_summary.py <chr1.hdf5/all_chr.hdf5> <pheno.filtered.hdf5> <correction_method> <peer.hdf5> <Kpop.hdf5> <covariates.hdf5> <window> <n_perm> <alpha> <nfolds> <samples_list> <summary_filename> '''


def getRelPos(pos,gene_pos):
	if gene_pos[3]==1:
		rv = pos-gene_pos[1]
	else:
		rv = gene_pos[2]-pos
	return rv

if __name__=='__main__':

	if len(sys.argv[1:]) < 11:
		sys.stderr.write('ERROR: missing parameters\n')
		usage()
		sys.exit(1)

	#load data from that are just needed to import data module 
	geno, pheno, correction_method, hdf_correction, Kpop, covariates,window,n_perm,alpha = sys.argv[1:10]
        window=float(window) #chromosomal window
	n_perm = int(n_perm) #number of permutations
	alpha = float(alpha) #level of significance for trans
	nfolds = int(sys.argv[10]) #number of folds used for gene expr matrix
	samples_list = sys.argv[11] #list of fold_j files
	summary = sys.argv[12] #outfile
	#populate dictionary with data for eqtl
	import eqtlsettings
	import data as DATA	
	#load doata
	data = DATA.data(geno, Kpop, pheno, covariates,  hdf_correction, correction_method,window)
	
	#check if file with samples list exist
	if os.path.isfile(samples_list) != True:
		sys.stderr.write('ERROR : '+samples_list+' not found\n')
		sys.exit(1)

	#load file with samples list
	fileslist = SP.loadtxt(samples_list,dtype='str') #TODO: these list of files needs to be produced before running this script in the make file

	#check if there are missing files according to nfolds		 
	if len(fileslist) != nfolds:
		sys.stderr.write('ERROR : missing files\n')
		sys.exit(1)

	#run the analysis					
	table = {}

	files = SP.sort(fileslist)

	#TODO: clean up the code here a little bit
	first = True
	for file in files:

		print file
		try:
			f = h5py.File(file,'r')
		except:
			print 'Oops! failed to open file'
                        sys.exit(1)
			continue

		geneIDs = f.keys()
		for geneID in geneIDs:
			fgene = f[geneID]
			try:
				temp = {}
				# store pv, pv_perm,cis qv,lambda,lambda_perm,lambda_empirical,beta
				#take the minimum pvalue for each gene. TODO: still don't know which strategy to use for trans!
				if window > 0: #is cis
					print "gene {0} kept".format(geneID)
					temp['geneID'] = SP.array([str(geneID)])
					temp['file'] = SP.array([str(file)])
					gene_pos = data.getGenePos(geneID)
					temp['gene_pos'] = gene_pos
					temp['lambda'] = fgene['lambda'][:,0]
                                	temp['lambda_perm'] = fgene['lambda_perm'][:,0]
					if n_perm > 1 :#if empirical pvalues have been computed			
						idx = fgene['pv_perm'][0,:].argmin()
						temp['pv_perm'] = fgene['pv_perm'][:,idx]
						temp['lambda_empirical']=fgene['lambda_empirical'][:,0]
					else:
						idx = fgene['qv'][0,:].argmin()
						temp['pv_perm'] = (SP.empty((1,))).astype(str)
						temp['pv_perm'][:] = 'NA' #no empirical pvalues has been computed
						temp['lambda_empirical']=fgene['lambda_empirical'][:]
					pos = fgene['pos'][[idx]]
					temp['gene_pos'] = gene_pos
					temp['chrom'] = fgene['chrom'][[idx]]
					temp['pos'] = pos
					temp['pv'] = fgene['pv'][:,idx]
					temp['qv'] = fgene['qv'][:,idx]
					temp['beta'] = fgene['beta'][:,idx]
				else: #is trans
					print "gene {0} kept".format(geneID)
					if n_perm > 1: #if empirical pvalues have been computed.
						idx = fgene['pv_perm'][0,:]<= alpha #boolean vector based on empirical pvalues
					else:
						idx = fgene['pv'][0,:]<= alpha #boolean vector based on nominal pvalues
					if sum(idx) > 0:
						temp['geneID'] = SP.tile(SP.array([str(geneID)]),sum(idx))
						temp['file'] = SP.tile(SP.array([str(file)]),sum(idx))
						gene_pos = data.getGenePos(geneID)
						temp['gene_pos'] = SP.tile(gene_pos,sum(idx))
						if n_perm > 1 : 
							temp['pv_perm'] = fgene['pv_perm'][0,:][idx]
							temp['lambda_empirical'] = SP.tile(fgene['lambda_empirical'][:,0],sum(idx))
						else:
							temp['pv_perm'] = (SP.empty((sum(idx),))).astype(str)
							temp['pv_perm'][:] = 'NA' #no empirical pvalues has been computed
							temp['lambda_empirical'] = SP.tile(fgene['lambda_empirical'][:],sum(idx))
						temp['pv'] = fgene['pv'][:][0,idx]
						temp['qv'] = fgene['qv'][:][0,idx]
						temp['beta'] = fgene['beta'][:][0,idx]
						pos = fgene['pos'][idx]
						temp['chrom'] = fgene['chrom'][idx]
						temp['pos'] = pos
						temp['lambda'] = SP.tile(fgene['lambda'][:,0],sum(idx))
						temp['lambda_perm'] = SP.tile(fgene['lambda_perm'][:,0],sum(idx))
					else:
						print "no pvalues <= alpha"
						pass
			except:
                                print "{0}:  nothing significant in here".format(geneID)
				pass
			#append the temp table into the big table
			for key in temp.keys():
				smartAppend(table,key,temp[key])
				#print table

		f.close()

	for key in table.keys():
		try: 
			table[key] = SP.concatenate(table[key])
		except:
			pass

        if len(table.keys())!=0:       
	# add corrected qvalues also across genes
                #table['qv_all'] = FDR.qvalues(table['qv'])
		pass
        else:
                print "Warning: no genes found"
	summary = h5py.File(summary,'w')
	smartDumpDictHdf5(table,summary)
	summary.close()

	sys.exit(0)

