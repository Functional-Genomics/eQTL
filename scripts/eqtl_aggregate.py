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
If cis, this script generates a final hdf5 file aggregating all the chr.hdf5 files and corrects the qvalues for multiple testing across all the genes tested; if trans it reads the unique summary file with trans results as input and corrects for multiple testing. 

Usage:

eqtl_aggregate.py <out.summary.cis.hdf5> <window> <n_perm>  <chr.lst>
'''
if __name__ == '__main__':

	if len(sys.argv[1:]) < 4:
		sys.stderr.write('ERROR: missing parameters\n')
		usage()
		sys.exit(1)

	outfile,window,n_perm = sys.argv[1:4]
	window = float(window)
	n_perm = int(n_perm)

	with open(sys.argv[4]) as f:
		chr = [i.strip() for i in f.readlines()]

	for file in chr:
		if os.path.isfile(file) != True:
			sys.stderr.write('ERROR: file '+file+' not found\n')
			sys.exit(1)

	if type(window) != (int) and  type(window) != (float):
		sys.stderr.write('\nERROR: Please use integer or float for window\n')
		usage()
		sys.exit(1)
	elif type(n_perm) != (int):
		sys.stderr.write('\nERROR: Please use integer for n_perm\n')
		usage()
		sys.exit(1)

	out = h5py.File(outfile,'w')
	#open a dictionary
	table = {}
	for file in chr:
		temp = {}
		i = h5py.File(file,'r')
		if i.keys() == []: #the file is empty
			print 'file {0} is empty\n'.format(i)
			pass
		else:
			temp['geneID'] =i['geneID'][:]
			temp['file'] = i['file'][:]
			temp['chrom'] = i['chrom'][:]
			temp['pos'] = i['pos'][:]
			temp['pv'] = i['pv'][:]
			temp['qv'] = i['qv'][:]
			temp['beta'] = i['beta'][:]
			temp['lambda'] = i['lambda'][:]
			temp['lambda_perm'] = i['lambda_perm'][:]
			temp['lambda_empirical'] = i['lambda_empirical'][:]
			temp['pv_perm'] = i['pv_perm'][:]				
		#append all file groups to the big table
		for key in temp.keys():
			smartAppend(table,key,temp[key])
	
	for key in table.keys():
		try:
			table[key] = sp.concatenate(table[key])
		except:
			pass

	#store the number of pvalues tested
	shape_pv = table['pv'].shape[0]
	#calculate qv_all,pv_perm_all
	if window != 0 and n_perm <=1: #if cis and no empirical pvalues
		table['qv_all'] = FDR.qvalues(table['qv'])
		table['window'] = [window] #add group with window size
		table['n_perm'] = [n_perm]
		table['pv_perm_all'] = sp.empty((shape_pv,))
		table['pv_perm_all'][:] = sp.NAN
	elif window !=0 and n_perm >1: # if cis and empirical pvalues
		table['qv_all'] = FDR.qvalues(table['qv'])
		table['pv_perm_all'] = FDR.qvalues(table['pv_perm'])
		table['window'] = [window]
		table['n_perm'] = [n_perm]
	elif window == 0 and n_perm <=1: #if trans and no empirical pvalues
		table['qv_all'] = FDR.qvalues(table['qv'])
		table['window'] = [window]
		table['n_perm'] = [n_perm]
		table['pv_perm_all'] = sp.empty((shape_pv,))
		table['pv_perm_all'][:] = sp.NAN
	else: #if trans and empirical pvalues
		table['qv_all'] = FDR.qvalues(table['qv'])
		table['pv_perm_all'] = FDR.qvalues(table['pv_perm'])
		table['window'] = [window]
		table['n_perm'] = [n_perm]

	#write within the file		
	smartDumpDictHdf5(table,out)
	out.close()
	sys.exit(0)

		


