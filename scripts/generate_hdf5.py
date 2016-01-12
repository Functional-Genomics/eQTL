#!/usr/bin/env python

import sys,os
import h5py
import scipy as sp
import pandas as pd

def usage():
	print '''
Usage:
generate_hdf.py <var_file.tsv> <geno_annotation.tsv> <chr> <chr.hdf5>
'''


if len(sys.argv[1:])<4:
	usage()
	sys.stderr.write('\nERROR: missing parameter\n')
	sys.exit(1)

file1,file2,chr,outfile=sys.argv[1:]

chr=str(chr)

if os.path.isfile(file1) != True:
	sys.stderr.write('\nERROR: file '+file1+' not found\n')
	sys.exit(1)
elif os.path.isfile(file2) != True:
	sys.stderr.write('\nERROR: file '+file2+' not found\n')
	sys.exit(1)


var_file=pd.read_csv(file1,sep='\t',index_col=0)
annotation = pd.read_csv(file2,sep='\t',index_col=0)

#bv=sp.in1d(annotation.index.values,var_file.index.values)

#check consistency between geno matrix and geno_annotation
#if sum(bv) != var_file.index.values.shape[0]:
#	bv2 = sp.in1d(var_file.index.values,annotation.index.values)
#	for var in var_file.index.values[~bv2]:
#		sys.stderr.write('\nERROR: variant '+var+' not found in the annotation file\n') 
#	sys.exit(1)

#extract variants annotations for a specific chr
chr_subset = annotation[annotation.ix[:,0]==chr]
chr_var_subset = chr_subset.index.values

var_file_subset = var_file.index.values[sp.in1d(var_file.index.values,chr_var_subset)]
#set array of geno values and transpose (samples X var)
matrix = var_file.values[sp.in1d(var_file.index.values,var_file_subset)].T.astype(float)
#set row_header
row_header= sp.array(var_file.columns.tolist())
#store indexes of the geno variant in the annotation list
i = map(lambda x:chr_var_subset.tolist().index(x),var_file_subset.tolist())
#set array of var positions and chromosome
pos = chr_subset.iloc[i,1].values.astype(float)
chrom = chr_subset.iloc[i,0].values.astype(int).astype(str)
#set array of allele. TODO: is empty for now.
allele = sp.zeros(pos.shape) 

#open an hdf5 file
hdf = h5py.File(outfile,'w')

#append the matrix, row_header,col_header (with subkeys)
dset = hdf.create_dataset('genotype/matrix',data=matrix)
dset = hdf.create_dataset('genotype/row_header/sample_ID', data=row_header)
dset = hdf.create_dataset('genotype/col_header/alleles',data=allele)
dset = hdf.create_dataset('genotype/col_header/chrom',data=chrom)
dset = hdf.create_dataset('genotype/col_header/pos',data=pos)

hdf.close()

sys.exit(0)


