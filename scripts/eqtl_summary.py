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
import pandas as pd

def usage():
	print '''
This script generates a final summary file from all the hdf5 files generated by the eqtl_cis/trans.py script. 
Usage:

eqtl_summary.py <chr1.hdf5/all_chr.hdf5> <pheno.filtered.hdf5> <correction_method> <peer.hdf5> <Kpop.hdf5> <covariates.hdf5> <window> <n_perm> <alpha> <nfolds> <samples_list> <summary_filename> <metainfo.tsv>'''


if __name__=='__main__':

	if len(sys.argv[1:]) < 12:
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
	summary = open(sys.argv[12],'a') #outfile
	metainfo = open(sys.argv[13],'w') #metainformation
	header = ['geneID','chrom','pos','pv','pv_perm','qv','beta','lambda','lambda_perm','file']
	header = pd.DataFrame(SP.array(header))
	header.T.to_csv(summary,mode = 'a',sep='\t',header=None,index=None)
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

	#store info about number of variants and phenotypeIDs and compute number of tests. Useful only for trans multiple test correction
	n_x = data.g['genotype/matrix'][:].shape[1] #number of variants
	n_y = data.p['phenotype/matrix'][:].shape[1] #number of phenotype_IDs
	n_tests = n_x*n_y #number of tests
	#run the analysis					

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
					temp['file'] = map(lambda x:x.split('/')[-1],SP.array([str(file)]))
					path = '\t'.join(SP.array([str(file)])[0].split('/')[:-1])
					temp['lambda'] = map(lambda x:round(x,3),fgene['lambda'][:,0])
                                	temp['lambda_perm'] = map(lambda x:round(x,3),fgene['lambda_perm'][:,0])
					if n_perm > 1 :#if empirical pvalues have been computed			
						idx = fgene['pv'][0,:].argmin()
						temp['pv_perm'] = map(lambda x:round(x,3),fgene['pv_perm'][:])
					else:
						idx = fgene['qv'][0,:].argmin()
						temp['pv_perm'] = (SP.empty((1,))).astype(str)
						temp['pv_perm'][:] = 'NA' #no empirical pvalues has been computed
					pos = fgene['pos'][[idx]]
					temp['chrom'] = fgene['chrom'][[idx]]
					temp['pos'] = pos
					temp['pv'] = map(lambda x:round(x,3),fgene['pv'][:,idx])
					try:
						temp['qv'] = map(lambda x:round(x,3),fgene['qv'][:,idx])
					except:
						temp['qv'] = fgene['qv'][:,idx]
					temp['beta'] = map(lambda x:round(x,3),fgene['beta'][:,idx])
				else: #is trans
					print "gene {0} kept".format(geneID)
					idx = fgene['pv'][0,:]<= alpha #boolean vector based on nominal pvalues
					s_idx = sum(idx)
					if s_idx > 0:
						temp['geneID'] = SP.tile(SP.array([str(geneID)]),s_idx)
						temp['file'] = map(lambda x:x.split('/')[-1],SP.tile(SP.array([str(file)]),s_idx))
						path = '\t'.join((SP.tile(SP.array([str(file)]),s_idx)[0]).split('/')[:-1])
						if n_perm > 1 : 
							temp['pv_perm'] = map(lambda x:round(x,3),SP.tile(fgene['pv_perm'][:,0],s_idx))
						else:
							temp['pv_perm'] = (SP.empty((s_idx,))).astype(str)
							temp['pv_perm'][:] = 'NA' #no empirical pvalues has been computed
						temp['pv'] = map(lambda x:round(x,3),fgene['pv'][:][0,idx])
						temp['qv'] = fgene['qv'][:][0,idx]
						temp['beta'] = map(lambda x:round(x,3),fgene['beta'][:][0,idx])
						pos = fgene['pos'][idx]
						temp['chrom'] = fgene['chrom'][idx]
						temp['pos'] = pos
						temp['lambda'] = map(lambda x:round(x,3),SP.tile(fgene['lambda'][:,0],s_idx))
						temp['lambda_perm'] = map(lambda x:round(x,3),SP.tile(fgene['lambda_perm'][:,0],s_idx))
					else:
						print "no pvalues <= alpha"
						pass
			except:
                                print "{0}:  nothing significant in here".format(geneID)
				pass
			temp_df = pd.DataFrame(SP.vstack((temp['geneID'][:],temp['chrom'][:],temp['pos'][:],temp['pv'][:],temp['pv_perm'][:],temp['qv'][:],temp['beta'][:],temp['lambda'][:],temp['lambda_perm'][:],temp['file'][:])).T)
			temp_df.to_csv(summary,mode = 'a',sep='\t',header=None,index=None)
		f.close()
		#write some meta information
	metainfo.write('n_tests\t'+str(n_tests)+'\n')
	metainfo.write('path\t'+path+'\n')
	metainfo.write('Number_of_variants_tested\t'+str(n_x)+'\n')
	metainfo.write('Number_of_phenotypesID_tested\t'+str(n_y)+'\n')
	metainfo.write('n_perm\t'+str(n_perm)+'\n')
	metainfo.write('pvalue_cutoff\t'+str(alpha)+'\n')
	metainfo.write('window\t'+str(window)+'\n')
	metainfo.write('correction_method\t'+correction_method+'\n')
		
	
	#close files
	metainfo.close()	
	summary.close()

	sys.exit(0)

