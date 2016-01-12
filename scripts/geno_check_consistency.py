#!/usr/bin/env python

import sys,os
import pandas as pd
import scipy as sp

def usage():
	print '''
Usage:
geno_check_consistency.py <var_file.tsv> <geno_annotation.tsv> 
'''


if len(sys.argv[1:])<2:
        usage()
        sys.stderr.write('\nERROR: missing parameter\n')
        sys.exit(1)

file1,file2=sys.argv[1:]

if os.path.isfile(file1) != True:
        sys.stderr.write('\nERROR: file '+file1+' not found\n')
        sys.exit(1)
elif os.path.isfile(file2) != True:
        sys.stderr.write('\nERROR: file '+file2+' not found\n')
        sys.exit(1)

#open the 2 files
var_file=pd.read_csv(file1,sep='\t',index_col=0)
annotation = pd.read_csv(file2,sep='\t',index_col=0)

#check if all the var in the var_file are in the geno_annotation file
bv=sp.in1d(var_file.index.values,annotation.index.values)

print 'Check consistency between the var_file and the annotation file'
#check consistency between var_file and geno_annotation file
if sum(bv) != var_file.index.values.shape[0]:
       for var in var_file.index.values[~bv]:
               sys.stderr.write('\nERROR: variant '+var+' not found in the annotation file\n')
       sys.exit(1)
else:
	print 'OK!'


sys.exit(0)

