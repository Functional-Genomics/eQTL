#!/usr/bin/env python


import h5py
import sys, os
import scipy as sp
import glob

def usage():
	print '''
This script generates a list of significant eqtls based on fdr threshold. The list contains the following fields:

1)GeneID
2)Chr
3)Position
4)p-value (nominal)
5)q-value (nominal, SNPs level correction)
6)q-value2 (nominal, gene level)
7)p-value_empirical (empirical, permutation) ['NA' if n_perm <=1 ]
8)q-value_empirical (empirical, gene level correction) ['NA' if n_perm <=1 ]
9)Effect size (nominal)
10)Lambda (nominal pvalues)
11)Lambda_perm (pvalues of the first permutation)


Usage:

get_results.py <fdr_threshold> <result.tsv> <final_res.hdf5>'''

def get_res(a,l,fdr):
	global genes,chromosome,position,pvalue,qvalue,qvalue_genes,pvalue_empirical,qvalue_empirical,beta,gen_control,gen_control_perm
	boolvector = (a<=fdr)
	genes = l[0][boolvector]
	chromosome = l[1][boolvector]
	chromosome = chromosome.astype(int)
	position = l[2][boolvector]
	position = position.astype(int)
	pvalue = l[3][boolvector]
	qvalue = l[4][boolvector]
	qvalue_genes = l[5][boolvector]
	pvalue_empirical=l[6][boolvector]
	qvalue_empirical=l[7][boolvector]
	beta =l[8][boolvector]
	gen_control=l[9][boolvector]
	gen_control_perm=l[10][boolvector]
	return genes,chromosome,position,pvalue,qvalue,qvalue_genes,pvalue_empirical,qvalue_empirical,beta,gen_control,gen_control_perm


if __name__ == '__main__':

	if len(sys.argv[1:]) < 3:
		sys.stderr.write('ERROR: missing parameters\n')
		usage()
		sys.exit(1)

	fdr,outfile = sys.argv[1:3]

	try:
		fdr = float(fdr)
	except:
		sys.stderr.write('ERROR: Please select integer of float value for fdr.\n')
		sys.exit(1)

	file = sys.argv[3]

	
	if os.path.isfile(file)!=True:
		sys.stderr.write('ERROR: file '+file+' not found\n')
		sys.exit(1)

	name_keys=['geneID','chrom','pos','pval','l_adj_pval','g_adj_pval','l_emp_pval','g_emp_adj_pval','beta','lambda_pval','lambda_perm','window','n_perm','file']
	#store header
	header="\t".join(name_keys[:11])

	#open tsv file
	out = open(outfile,'w')
	out.write(header+'\n')

	#open hdf5 file with final results
	i_file = h5py.File(file,'r')

	#check if file is non-empty
	if i_file.keys() == [] : #file is empty:
		out.close()
		sys.stdout.write('ERROR: file '+file+' is empty\n')
		sys.exit(1)
	
	#check if file has all the expected keys
	for key in name_keys:
		if key not in i_file: 
			std.error.write('ERROR: key '+key+' is not present in file '+i_file+'\n')
			sys.exit(1)

	#if everything passed numpy array check
	n_perm = i_file['n_perm'][0] #n_perm
	if n_perm <= 1: #no empirical pvalues
		a = i_file['g_adj_pval'][:] #use qv_all to select genes based on fdr
	else:
		a = i_file['g_emp_adj_pval'][:] #use pv_perm_all to select genes based on fdr

	for x in xrange(len(name_keys)):
		name_keys[x] = i_file[name_keys[x]][:] #assign matrix to each variable

	o = get_res(a,name_keys,fdr)
	if genes.shape[0] == 0:
		print 'no significant result found for file {0}'.format(outfile)
		pass
	else:	
		#initialise an array with final results
		finalarray = sp.arange(o[0].shape[0])
		finalarray = o[0]
		for res in o[1:]:
			finalarray = sp.column_stack((finalarray,res))
		#write into the outfile
		for line in finalarray:
			r = '\t'.join(line)
			out.write(r+'\n')

	out.close()
	sys.exit(0)

		


