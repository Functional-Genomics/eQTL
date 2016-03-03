#!/usr/bin/env python

import sys,os
import h5py
import scipy as sp
from utils import corr_matrix
import pandas as pd

def usage():
	print '''
This scripts substitutes the mutation burden in 0/1 and generates an hdf5 file using the matrix of genotype veriants and the annotation bed file.
 
Usage:
cat <var_filtered_stdin> | generate_hdf.py <var_annotation.bed> <chr> <chr.hdf5> 
Optional parameter: <skip_kinship>
'''


if len(sys.argv[1:])<3:
	usage()
	sys.stderr.write('\nERROR: missing parameter\n')
	sys.exit(1)

if len(sys.argv[1:]) == 3:
	file2,chr,outfile=sys.argv[1:]
	skip_kinship=False #generate kinship by default
elif len(sys.argv[1:]) == 4:
	file2,chr,outfile,skip_kinship = sys.argv[1:]
	skip_kinship=True #generate kinship with NaN

chr=str(chr)

if os.path.isfile(file2) != True:
	sys.stderr.write('\nERROR: file '+file2+' not found\n')
	sys.exit(1)

#open an hdf5 file
hdf = h5py.File(outfile,'w')

sys.stderr.write('\nReading from stdin... ')
file1 = sys.stdin
var_file=pd.read_csv(file1,sep='\t',index_col=0)
sys.stderr.write('complete.\n')

sys.stderr.write('Reading from '+file2+'... ')
annotation = pd.read_csv(file2,sep='\t',index_col=0,header=None)
sys.stderr.write('done.')

chr_subset = annotation[annotation.index.values==chr]
chr_var_subset = chr_subset[3].values #take names from bed file

var_file_subset = var_file.index.values[sp.in1d(var_file.index.values,chr_var_subset)]
#set array of geno values and transpose (samples X var)
matrix = var_file.values[sp.in1d(var_file.index.values,var_file_subset)].T.astype(float)
#save a copy of the matrix with burden
dset = hdf.create_dataset('genotype/burden_matrix',data=matrix)
burden = matrix.copy()
#compute kinship based on optional parameter (skip_kinship)
matrix,K = corr_matrix(matrix[:],skip_kinship=skip_kinship) #this step overwrite the matrix array (standardised genotypes)
#change mutation burden into 1
burden[burden>=1]=1
#set row_header
row_header= sp.array(var_file.columns.tolist())
#store indexes of the geno variant in the annotation list
i = map(lambda x:chr_var_subset.tolist().index(x),var_file_subset.tolist())
#set array of var positions (mean point of the annotation interval) and chromosome
pos = sp.mean(chr_subset.iloc[i,[0,1]],axis=1).values.astype(float)
chrom = chr_subset.index[i].values.astype(int).astype(str)
#set array of allele. TODO: is empty for now.
allele = sp.zeros(pos.shape) 

#append the matrix, row_header,col_header (with subkeys)
dset = hdf.create_dataset('genotype/matrix',data=burden)
dset = hdf.create_dataset('genotype/Kpop',data=K)
dset = hdf.create_dataset('genotype/row_header/sample_ID', data=row_header)
dset = hdf.create_dataset('genotype/col_header/alleles',data=allele)
dset = hdf.create_dataset('genotype/col_header/chrom',data=chrom)
dset = hdf.create_dataset('genotype/col_header/pos',data=pos)

hdf.close()

sys.exit(0)


