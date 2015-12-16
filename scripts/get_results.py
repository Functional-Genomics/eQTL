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

def get_res(a,b,c,d,e,f,fdr,l,lp):
	global genes,chromosome,position,qvalue,qvalue_genes,beta,gen_control,gen_control_perm
	boolvector = (e<=fdr)
	genes = a[boolvector]
	chromosome = b[boolvector]
	chromosome = chromosome.astype(int)
	position = c[boolvector]
	position = position.astype(int)
	qvalue = d[boolvector]
	qvalue_genes = e[boolvector]
	beta =f[boolvector]
	gen_control=l[boolvector]
	gen_control_perm=lp[boolvector]
	return genes,chromosome,position,qvalue,qvalue_genes,beta,gen_control,gen_control_perm

def check_numpy_df(hdf5file,key):
	try:
		hdf5file[key][:]
		o=0
	except:
		o=1
	return o


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

	#open tsv file
	out = open(outfile,'w')

	#open hdf5 file with final results
	i = h5py.File(file,'r')

	#check if file is non-empty
	if i.keys() == [] : #file is empty:
		std.error.write('ERROR: file '+file+' is empty\n')
		sys.exit(1)
	
	name_keys=['geneID','chrom','pos','qv','qv_all','beta','lambda','lambda_perm','window']
	
	#check if file has all the expected keys
	for key in name_keys:
		o=check_numpy_df(i,key)
		if o == 1:
			std.error.write('ERROR: key '+key+' is not present in file '+file+'\n')
			sys.exit(1)

	#if everything passed numpy array check

	a = i[name_keys[0]][:] #geneID
	b = i[name_keys[1]][:] #chrom
	c = i[name_keys[2]][:] #pos
	d = i[name_keys[3]][:] #qv
	e = i[name_keys[4]][:] #qv_all
	f = i[name_keys[5]][:] #beta
	window = i[name_keys[-1]][0] #window

	if window == 0: # is trans
		l=i[name_keys[6]][:][0] # TODO: check why this applies only to trans
		lp=i[name_keys[7]][:][0] 
	else: #is cis
		l=i[name_keys[6]][:]
		lp=i[name_keys[7]][:]

	o = get_res(a,b,c,d,e,f,fdr,l,lp)
	if genes.shape[0] == 0:
		print 'no significant result found for file {0}'.format(outfile)
		pass
	else:	
		#initialise an array with final results
		finalarray = sp.arange(o[0].shape[0])
		finalarray = o[0]
		for res in o[1:]:
			finalarray = sp.column_stack((finalarray,res))
		out.write('\t'.join(name_keys[:-1])+'\n')
		for line in finalarray:
			r = '\t'.join(line)
			out.write(r+'\n')

	out.close()
	sys.exit(0)

		


