#!/usr/bin/env python


import h5py
import sys, os
import scipy as sp
import glob
import pandas as pd
import gzip

def usage():
	print '''
This script generates a list of significant eqtls based on fdr threshold. The list contains the following fields:

1)GeneID
2)Chr
3)Position
4)Variant name
5)p-value (nominal)
6)q-value (nominal, SNPs level correction)
7)q-value2 (nominal, gene level)
8)p-value_empirical (empirical, permutation) ['NA' if n_perm <=1 ]
9)q-value_empirical (empirical, gene level correction) ['NA' if n_perm <=1 ]
10)Effect size (nominal)
11)Lambda (genomic control based on nominal pvalues)
12)Lambda_perm_mean (mean value of genomic control based on n permutation)
13)Lambda_perm_std (std value of genomic control based on n permutation)


Usage:

eqtl_step4.py <file_res.tsv> <metainfo.tsv> <fdr_threshold> <outfile.tsv>'''

def get_res(a,l,fdr):
	global genes,chromosome,position,var_name,pvalue,qvalue,qvalue_genes,pvalue_empirical,qvalue_empirical,beta,gen_control,gen_control_perm_mean,gen_control_perm_std
	boolvector = (a<=fdr)
	genes = l[0][boolvector]
	chromosome = l[1][boolvector]
	chromosome = chromosome.astype(int)
	position = l[2][boolvector]
	position = position.astype(int)
	var_name = l[3][boolvector]
	pvalue = l[4][boolvector]
	qvalue = l[5][boolvector]
	qvalue_genes = l[6][boolvector]
	pvalue_empirical=l[7][boolvector]
	qvalue_empirical=l[8][boolvector]
	beta =l[9][boolvector]
	gen_control=l[10][boolvector]
	gen_control_perm_mean=l[11][boolvector]
	gen_control_perm_std = l[12][boolvector]
	return genes,chromosome,position,var_name,pvalue,qvalue,qvalue_genes,pvalue_empirical,qvalue_empirical,beta,gen_control,gen_control_perm_mean,gen_control_perm_std


if __name__ == '__main__':

	if len(sys.argv[1:]) < 4:
		sys.stderr.write('ERROR: missing parameters\n')
		usage()
		sys.exit(1)

	fdr,outfile = sys.argv[3:5]

	try:
		fdr = float(fdr)
	except:
		sys.stderr.write('ERROR: Please select integer of float value for fdr.\n')
		sys.exit(1)

	file = sys.argv[1]
	metainfo = sys.argv[2]
	
	if os.path.isfile(file)!=True:
		sys.stderr.write('ERROR: file '+file+' not found\n')
		sys.exit(1)
	if os.path.isfile(metainfo)!=True:
		sys.stderr.write('ERROR: file '+metainfo+' not found\n')
		sys.exit(1)

	name_keys=['geneID','chrom','pos','var_name','pval','l_adj_pval','g_adj_pval','l_emp_pval','g_emp_adj_pval','beta','lambda_pval','lambda_perm_mean','lambda_perm_std']
	header = "\t".join(name_keys)+'\n'
	#open tsv file
	out = gzip.open(outfile,'wb')
	metainfo = pd.read_csv(metainfo,sep="\t",header=None,index_col=[0])

	#open hdf5 file with final results
#	i_file = h5py.File(file,'r')
	i_file = pd.read_csv(file,sep="\t",na_values='NA',compression='gzip')
	#check if file is non-empty
	if i_file.shape[0] == 1 : #file is empty:
		out.close()
		sys.stdout.write('WARNING: file '+file+' is empty\n')
		sys.exit(0)
	
	#check if file has all the expected keys
	for key in name_keys:
		if key not in i_file.columns:
			sys.stderr.write('ERROR: key '+key+' is not present in file '+i_file+'\n')
			sys.exit(1)
	
	#if everything passed numpy array check
	n_perm = int(metainfo[metainfo.index=='n_perm'].values[0][0]) #n_perm
	window = metainfo[metainfo.index=='window'].values[0][0]
	if n_perm <= 100: #no empirical pvalues
		a = i_file['g_adj_pval'].values #use qv_all to select genes based on fdr
	else:
		a = i_file['g_emp_adj_pval'].values #use pv_perm_all to select genes based on fdr
	for x in xrange(len(name_keys)):
		name_keys[x] = i_file[name_keys[x]].values #assign matrix to each variable

	o = get_res(a,name_keys,fdr)
	if genes.shape[0] == 0:
		print 'no significant result found for file {0}'.format(outfile)
		out.write(header)
		pass
	else:	

		out.write(header)
		#initialise an array with final results
		finalarray = sp.arange(o[0].shape[0])
		finalarray = o[0]
		for res in o[1:]:
			finalarray = sp.column_stack((finalarray,res))
		#write into the outfile
		for line in finalarray:
			line[line.astype(str)=='nan']='NA'
			line=line.astype(str)
			r = '\t'.join(line)
			out.write(r+'\n')

	out.close()
	sys.exit(0)

		


