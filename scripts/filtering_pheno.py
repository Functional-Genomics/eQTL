#!/usr/bin/env python

import sys, os
import h5py
import scipy as SP
import pandas as pd
import normalise as NM


def usage():
	print """ 
This script filters the matrix of gene expression (samples X genes) values based on a gene expression threshold (e.g. FPKMs) and minimum percentage n of samples meeting this threshold. It outputs a tsv file with (genes X samples)

Usage:
filtering_pheno.py <pheno.matched.tsv> <min_expr> <min_perc_samples> <pheno.filtered.tsv> <transpose>


- min_exp = [INT | FLOAT]
- min_perc_samples = [INT | FLOAT]
- transpose [OPTIONAL] = [y | n; DEFAULT: y]
"""

#check arguments
if len(sys.argv[1:]) < 4 :
	sys.stderr.write("ERROR: missing parameter\n")
	usage()
	sys.exit(1)
elif len(sys.argv[1:]) == 5:
	pheno,threshold,perc_samples,phenout,transpose=sys.argv[1:]
else:
	pheno,threshold,perc_samples,phenout=sys.argv[1:]
	transpose = 'y'

#check files/arguments exist
if os.path.isfile(pheno) != True:
	sys.stderr.write("ERROR: file "+pheno+" not found\n")
	usage()
	sys.exit(1)
#elif expr_transform not in [ 'gauss', 'log', 'none' ]:
#        sys.stderr.write("ERROR: method "+expr_transform+" not available\n")
#        usage()
#        sys.exit(1)


pheno = pd.read_csv(pheno, sep='\t', index_col=[0])
threshold = float(threshold)
perc_samples = float(perc_samples)

#get number of genes and samples and matrix of gene expression
samples = pheno.shape[0]
#apply filters
samples_threshold=(perc_samples/100.0)*samples
vector=SP.sum((pheno.values[:]>=threshold), axis=0)
booleanvector=(vector>=samples_threshold)
indexes=SP.where(booleanvector == True)[0]

print "{0} genes retained after filtering".format(len(indexes))

#filter original matrix of gene expression
filtered_matrix = pheno.iloc[:,indexes]

#output the matrix 
if transpose == 'y':
	filtered_matrix.T.to_csv(phenout, sep='\t', index=True)
else:
	filtered_matrix.to_csv(phenout, sep='\t', index=True)

#apply transformation to gene expression data
#if expr_transform == 'gauss':
#	Y=NM.gaussianize(filtered_matrix[:])
#elif expr_transform == 'log':
#	Y=NM.logtransform(filtered_matrix[:])
#else:
#	Y=filtered_matrix[:]


sys.exit(0)
