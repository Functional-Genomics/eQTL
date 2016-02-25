#!/usr/bin/env python

import sys, os
import h5py
import scipy as SP
import pandas as pd
import normalise as NM


def usage():
	print """ 
This script optionally transforms filtered gene counts to gauss distribution. It assumes that the matrix is (samples X genes). If not please set transpose = y.

Usage:
filtering_pheno.py <pheno.filtered.qn.tsv> <expr_transform> <pheno.filtered.qn.normal.tsv> <transpose>

- expr_transform = [gauss | log | none]
- transpose = [y | n]
"""

#check arguments
if len(sys.argv[1:]) < 4:
	sys.stderr.write("ERROR: missing parameter\n")
	usage()
	sys.exit(1)

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

if transpose == 'y':
	pheno = pheno.T[:]
	


#apply transformation to gene expression data
if expr_transform == 'gauss':
	Y=NM.gaussianize(pheno.values[:])
elif expr_transform == 'log':
	Y=NM.logtransform(pheno.values[:])
else:
	Y=pheno.values[:]

#substitute values
pheno.values[:] = Y[:]
#write to output
pheno.to_csv(phenout,sep='\t',index=True)

sys.exit(0)
