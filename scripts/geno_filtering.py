#!/usr/bin/env python 

import sys,os
import pandas as pd
import scipy as sp

def usage():
	print ''' 
This script selects genomic variants based on the selected fraction of samples carrying the variant

Usage: cat matrix.tsv | geno_filtering.py <sample_fraction_threshold>
'''

if len(sys.argv[1:]) < 1:
	usage()
	sys.stderr.write('\nERROR: missing parameter\n')
	sys.exit(1)

threshold = sys.argv[1]

try:
	threshold = float(threshold)
except:
	usage()
	sys.stderr.write('\nERROR: threshold must be float\n')
	sys.exit(1)

#read file from stdin
file1 = sys.stdin
#read the tsv file
matrix = pd.read_csv(file1,sep='\t',index_col=0)
#change index type into string 
if matrix.index.dtype == (float):
	matrix.index=matrix.index.values.astype(int).astype(str)
elif matrix.index.dtype == (int):
	matrix.index=matrix.index.values.astype(str)

#select rows based on number of samples with at least 1 event (it assumes a SNPs by Samples matrix)
samples = matrix.shape[1]
sample_threshold = threshold * samples
filt_matrix = matrix[matrix.gt(0,axis='rows').sum(1) >= sample_threshold]

#write to stdout
filt_matrix.to_csv(sys.stdout,sep='\t',header=True,index=True)

sys.exit(0)

