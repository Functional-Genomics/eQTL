#!/usr/bin/env python


import h5py
import sys, os
import scipy as sp
import glob
from utils import smartAppend
import limix.stats.fdr as FDR
from utils import smartDumpDictHdf5
from utils import getLambda
import readline # because of an error with importing rpy2.robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import pandas as pd


stats = importr('stats')

def usage():
	print '''
If cis, this script generates a final hdf5 file aggregating all the chr.hdf5 files and corrects the qvalues for multiple testing across all the genes tested; if trans it reads the unique summary file with trans results as input and corrects for multiple testing using Bonferroni correction. 

Usage:

eqtl_aggregate.py <summary.fold.tsv> <metainfo.tsv> <outfile.tsv> 
'''
if __name__ == '__main__':

	if len(sys.argv[1:]) < 3:
		sys.stderr.write('ERROR: missing parameters\n')
		usage()
		sys.exit(1)

	file,metainfo,outfile = sys.argv[1:]
#	with open(file) as f:
#		chr = [i.strip() for i in f.readlines()]	
	if os.path.isfile(file) != True:
		sys.stderr.write('ERROR: file '+file+' not found\n')
		sys.exit(1)
	elif os.path.isfile(metainfo)!= True:
		sys.stderr.write('ERROR: file '+metainfo+' not found\n')
		sys.exit(1)

	file = pd.read_csv(file,sep='\t')
	metainfo = pd.read_csv(metainfo,sep='\t',index_col=[0],header=None)

	n_tests = int(metainfo[metainfo.index == 'n_tests'].values[0][0])
	window = float(metainfo[metainfo.index == 'window'].values[0][0])
	n_genes = int(metainfo[metainfo.index == 'Number_of_phenotypesID_tested'].values[0][0])
	n_perm = int(metainfo[metainfo.index == 'n_perm'].values[0][0])
	#open out tsv file
	out = open(outfile,'w')
	if file.shape[0] == 1:
		sys.stdout.write('WARNING: file {0} is empty\n'.format(i))
	else:
		pval = file['pv'].astype(float)
#		l_adj_pval = file['qv'].astype(float)
#		beta = file['beta'].astype(float)
#		lambda_pval= file['lambda'].astype(float)
#		lambda_perm= file['lambda_perm'].astype(float)
#		l_emp_pval = file['pv_perm'].astype(float)
		#store the number of pvalues tested
		shape_pv = pval.shape[0]
		#calculate qv_all,pv_perm_all
		if window != 0 and n_perm <=1: #if cis and no empirical pvalues
			g_adj_pval = FDR.qvalues(file['qv'].astype(float),m=n_genes)
			g_emp_adj_pval = (sp.empty((shape_pv,))).astype(str)
			g_emp_adj_pval[:] = 'NA' #fill an empty array with NA values
		elif window !=0 and n_perm >1: # if cis and empirical pvalues
			g_adj_pval = FDR.qvalues(file['qv'].astype(float),m=n_genes)
			g_emp_adj_pval = FDR.qvalues(file['pv_perm'].astype(float),m=n_genes)
		elif window == 0 and n_perm <=1: #if trans and no empirical pvalues
			pval = file['pv'].astype(float)
			g_adj_pval = sp.array(stats.p_adjust(FloatVector(pval.tolist()),method = 'bonferroni',n=float(n_tests))) #compute bonferroni adjusted across nominal pvalues
			g_emp_adj_pval = (sp.empty((shape_pv,))).astype(str) #fill an empty array with NA values
			g_emp_adj_pval[:] = 'NA'
		else: #if trans and empirical pvalues
			g_adj_pval = (sp.empty((shape_pv,))).astype(str) #fill an empty array with NA values
			g_emp_adj_pval= sp.array(stats.p_adjust(FloatVector((file['pv_perm'].astype(float)).tolist()),method = 'bonferroni',n=float(n_tests))) #compute bonferroni adjusted across empirical pvalues
	file.insert(0,'g_adj_pval',pd.DataFrame(g_adj_pval))
	file.insert(0,'g_emp_adj_pval',pd.DataFrame(g_emp_adj_pval))
	file = file[['geneID','chrom','pos','pv','pv_perm','g_emp_adj_pval','qv','g_adj_pval',	'beta','lambda','lambda_perm','file']]
	file.rename(columns={'pv':'pval','pv_perm':'l_emp_pval','qv':'l_adj_pval','lambda':'lambda_pval'},inplace=True)
	file.to_csv(out,sep='\t',header=True,index=None,na_rep='NA')
	out.close()	
	

	sys.exit(0)

		


