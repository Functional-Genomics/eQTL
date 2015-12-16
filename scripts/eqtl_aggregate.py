#!/usr/bin/env python


import h5py
import sys, os
import scipy as sp
import glob
from utils import smartAppend
import limix.stats.fdr as FDR
from utils import smartDumpDictHdf5
from utils import getLambda

def usage():
	print '''
This script generates a final hdf5 file aggregating all the chr.hdf5 files. It corrects the qvalues for multiple testing across all the genes tested.

Usage:

eqtl_aggregate.py <summary.hdf5> <window>  <1.hdf5> [<2.hdf5> ... ]
'''
if __name__ == '__main__':

	if len(sys.argv[1:]) < 3:
		sys.stderr.write('ERROR: missing parameters\n')
		usage()
		sys.exit(1)

	outfile,window = sys.argv[1:3]
	window = float(window)

	chr = sys.argv[3:]

	for file in chr:
		if os.path.isfile(file)!=True:
			sys.stderr.write('ERROR: file '+file+' not found\n')
			sys.exit(1)

	out = h5py.File(outfile,'w')
	#open a dictionary
	table = {}
	for file in chr:
		i = h5py.File(file,'r')
		if i.keys() == []: #the file is empty
			print 'file {0} is empty\n'.format(i)
			pass
		else:
			if window != 0 : # is cis:
				temp = {}
				temp['geneID'] =i['geneID'][:]
				temp['chrom'] = i['chrom'][:]
				temp['pos'] = i['pos'][:]
				temp['pv'] = i['pv'][:]
				temp['qv'] = i['qv'][:]
				temp['beta'] = i['beta'][:]
				temp['lambda'] = i['lambda'][:]
				temp['lambda_perm'] = i['lambda_perm'][:]
			else: #is trans
				temp = {}
				table['geneID'] =i['geneID'][:]
				temp['chrom'] = i['chrom'][:]
				temp['pos'] = i['pos'][:]
				temp['pv'] = i['pv'][:]
				temp['beta'] = i['beta'][:]
				temp['pv_perm'] = i['pv_perm'][:]
		#append all file groups to the big table
		for key in temp.keys():
			smartAppend(table,key,temp[key])
	
	for key in table.keys():
		try:
			table[key] = sp.concatenate(table[key])
		except:
			pass
	
	#calculate qv_all
	if window != 0:
		table['qv_all'] = FDR.qvalues(table['qv'])
	else: #if it's trans
		pvalues=sp.hstack(table['pv'][0:])
		beta = sp.hstack(table['beta'][0:])
		pvalues_perm=sp.hstack(table['pv_perm'][0:])
		pos = table['pos'][:]
		chrom = table['chrom'][:]
		del table['pv']
		del table['pos']
		del table['chrom']
		del table['beta']
		del table['pv_perm']
		l=getLambda(pvalues)[0]
		lp=getLambda(pvalues_perm)[0]
		smartAppend(table,'lambda',l)
		smartAppend(table,'lambda_perm',lp)
		for gene_pv in xrange(len(pvalues)):
			#caclulate qvalues correcting across all the snps tested per gene
			qv=FDR.qvalues(pvalues[gene_pv]) 
			idx= qv.argmin() #tale the index of the minumum pvalue per gene
			#append to the big table of trans eqtls
			smartAppend(table,'qv',qv[idx]) #append minimum qvalues 
			smartAppend(table,'pv',pvalues[gene_pv][idx]) #append minimum pvalue per gene
			smartAppend(table,'pos',pos[idx]) #append the position of the SNP with the min qvalue per gene
			smartAppend(table,'chrom',chrom[idx]) #append the chrom of the SNP with the min qvalue per gene
			smartAppend(table,'beta',beta[gene_pv][idx]) #append the beta for the SNP with the min qvalue per gene
		table['qv_all'] = FDR.qvalues(sp.array(table['qv'])) #if you use qvalues method does not work. TODO: understand what qvalues1 is doing

	#write within the file		
	smartDumpDictHdf5(table,out)
	out.close()
	sys.exit(0)

		


