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
3)Position(Leading Variant)
4)p-value (Leading variant)
5)q-value (adjusted p-value of the leading variant)
5)Effect size
6)Standard error for all variants #???


Usage:

get_results.py <fdr_threshold> <result.tsv> <final_res.hdf5>'''

def get_res(a,b,c,d,e,f,g,fdr,l,lp):
	global genes,chromosome,position,pvalue,qvalue,beta,gen_control,gen_control_perm
	boolvector = (e<=fdr)
	genes = a[boolvector]
	chromosome = b[boolvector]
	chromosome = chromosome.astype(int)
	position = c[boolvector]
	position = position.astype(int)
	pvalue = d[boolvector]
	qvalue_local = e[boolvector]
	qvalue_genes = f[boolvector]
	beta =g[boolvector]
	gen_control=l[boolvector]
	gen_control_perm=lp[boolvector]
	return genes,chromosome,position,pvalue,qvalue_local,qvalue_genes,beta,gen_control,gen_control_perm

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

	out = open(outfile,'w')
	out.write('geneID\tchrom\tpos\tpv\tqv\tqv_all\tbeta\tlambda\tlambda_perm\n')


	i = h5py.File(file,'r')
	a = i['geneID'][:]
	b = i['chrom'][:]
	c = i['pos'][:]
	d = i['pv'][:]
	e = i['qv'][:]
	f = i['qv_all'][:]
	g = i['beta'][:]
	window = i['window'][0]

	if window == 0: # is trans
		l=i['lambda'][:][0] # TODO: check why this applies only to trans
		lp=i['lambda_perm'][:][0] 
	else: #is cis
		l=i['lambda'][:]
		lp=i['lambda_perm'][:]

	o = get_res(a,b,c,d,e,f,g,fdr,l,lp)
	if genes.shape[0] == 0:
		print 'no significant result found for file {0}'.format(outfile)
		pass
	else:	
		#initialise an array with final results
		finalarray = sp.arange(o[0].shape[0])
		finalarray = o[0]
		for res in o[1:]:
			finalarray = sp.column_stack((finalarray,res))
		for line in finalarray:
			r = '\t'.join(line)
			out.write(r+'\n')

	out.close()
	sys.exit(0)

		


