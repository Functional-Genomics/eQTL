#!/usr/bin/env python


import h5py
import sys,os
import scipy as sp


def usage():
	print ''' 
\nThis script aggregates different chromosome matrices into the same hdf5 

Usage:
aggregate_chromosomes.py <chr_list.lst> <outfile.hdf5>
'''

def stack_array(m,n):
	''' add column to a numpy array'''
	m = sp.hstack((m,n))
	return m

if __name__ == '__main__':
	if len(sys.argv[1:])<2:
		sys.stderr.write('ERROR: missing parameters\n')
		usage()
		sys.exit(1)

	#read arguments
	chr_list,outfile = sys.argv[1:]
	#open list
	chr_list = open(chr_list,'r')
	chr_list = chr_list.readlines()
	chr_list = [i.strip() for i in chr_list]
	#initialise some empty arrays
	matrix = ''
	chromosome = ''
	pos = ''
	var_name = ''
	samples = ''
	alleles = ''

	#check is some arguments are missed
	for chr in chr_list:
		if os.path.isfile(chr) != True:
			sys.stderr.write('ERROR: file '+chr+' not found\n')
			sys.exit(1)
		else:
			chr = h5py.File(chr,'r')
			if matrix == '' and chromosome == '' and pos == '' and var_name == '' and samples == '' and alleles == '' :
				matrix,chromosome,pos,var_name,samples,alleles = chr['genotype/matrix'][:],chr['genotype/col_header/chrom'][:],chr['genotype/col_header/pos'][:],chr['genotype/col_header/var_names'][:],chr['genotype/row_header/sample_ID'][:],chr['genotype/col_header/alleles'][:]
			else:
				matrix = stack_array(matrix,chr['genotype/matrix'][:])
				chromosome = stack_array(chromosome,chr['genotype/col_header/chrom'][:])
				pos = stack_array(pos,chr['genotype/col_header/pos'][:])
				var_name = stack_array(var_name,chr['genotype/col_header/var_names'][:])
				if alleles.shape[0] == 0:
					sys.stderr.write('No significant genomic variants found for this chromosome!\n')
					pass
				else:
					alleles = sp.concatenate((alleles,chr['genotype/col_header/alleles'][:]))
					


	outfile = h5py.File(outfile,'w')
	dset = outfile.create_dataset('genotype/matrix',data = matrix[:])
	dset = outfile.create_dataset('genotype/row_header/sample_ID',data = samples[:])
	dset = outfile.create_dataset('genotype/col_header/pos',data = pos[:])
	dset = outfile.create_dataset('genotype/col_header/chrom',data = chromosome[:])
	dset = outfile.create_dataset('genotype/col_header/alleles',data = alleles[:])
	dset = outfile.create_dataset('genotype/col_header/var_names',data = var_name[:])

	outfile.close()
	sys.exit(0)

