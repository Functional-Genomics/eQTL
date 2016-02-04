#!/usr/bin/env python 

import sys,os
import pandas as pd
import scipy as sp

def usage():
	print '''

Usage:
geno_filtering.py <stdin> <threshold> <filtered_matrix_stdout>
'''

if len(sys.argv[1:]) < 1:
	usage()
	sys.stderr.write('\nERROR: missing parameter\n')
	sys.exit(1)

threshold = sys.argv[1]


#if os.path.isfile(file1) != True:
#	sys.stderr.write('\nERROR: file '+file1+' not found\n')
#	sys.exit(1)

try:
	threshold = int(threshold)
except:
	usage()
	sys.stderr.write('\nERROR: threshold must be int\n')
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

#select rows based on number of samples with at least 1 event
filt_matrix = matrix[matrix.gt(0,axis='rows').sum(1) >= threshold]

filt_matrix.to_csv(sys.stdout,sep='\t',header=True,index=True)

sys.exit(0)

