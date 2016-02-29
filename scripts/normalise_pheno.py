#!/usr/bin/env python

import sys, os
import h5py
import scipy as SP
import pandas as pd
import normalise as NM


def usage():
	print """ 
This script optionally transforms filtered gene counts to gauss distribution. It takes a matrix of genes X samples and transposes the final matrix.

Usage:
normalising_pheno.py <pheno.filtered.qn.tsv> <expr_transform> <pheno.filtered.qn.normal.tsv> <transpose>

- expr_transform = [gauss | log | none]
- OPTIONAL : transpose = [y | n] 
"""

#check arguments
if len(sys.argv[1:]) < 3:
	sys.stderr.write("ERROR: missing parameter\n")
	usage()
	sys.exit(1)

if len(sys.argv[1:]) == 3:
	pheno,expr_transform,phenout=sys.argv[1:]
	transpose = 'y'
elif len(sys.argv[1:]) == 4:
	pheno,expr_transform,phenout,transpose=sys.argv[1:]
	
#check files/arguments exist
if os.path.isfile(pheno) != True:
	sys.stderr.write("ERROR: file "+pheno+" not found\n")
	usage()
	sys.exit(1)
elif expr_transform not in [ 'gauss', 'log', 'none' ]:
        sys.stderr.write("ERROR: method "+expr_transform+" not available\n")
        usage()
        sys.exit(1)
elif transpose != 'n' and transpose != 'y':
	sys.stderr.write("ERROR: please use \'y\' or \'n\' for transpose option\n")
	usage()
	sys.exit(1)


pheno = pd.read_csv(pheno, sep='\t', index_col=[0])
	
#apply transformation to transposed gene expression data
if expr_transform == 'gauss':
	Y=NM.gaussianize(pheno.T.values[:])
elif expr_transform == 'log':
	Y=NM.logtransform(pheno.T.values[:])
else:
	Y=pheno.values[:]


#substitute values
pheno = pheno.T[:]
pheno.values[:] = Y[:]

if transpose == 'n':
	pheno = pheno.T[:] #to have a matrix of genes X samples
else:
	pass #to have a matrix of samples X genes

#write to output
pheno.to_csv(phenout,sep='\t',index=True)

sys.exit(0)
