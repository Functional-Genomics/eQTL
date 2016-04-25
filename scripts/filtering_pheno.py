#!/usr/bin/env python

import sys, os
import h5py
import scipy as SP
import pandas as pd
from utils import select_genes_on_SNR as SNR


def usage():
	print """ 
This script filters the matrix of gene expression (genes X samples) values based on a gene expression threshold (e.g. FPKMs) and minimum percentage n of samples meeting this threshold. It outputs a tsv file with (genes X samples)

Usage:
filtering_pheno.py <pheno.matched.tsv> <min_expr> <min_perc_samples> <snr_threshold> <pheno.filtered.tsv>


- min_exp = [INT | FLOAT]
- min_perc_samples = [INT | FLOAT]
- snr_threshold = [INT]
"""

#check arguments
if len(sys.argv[1:]) < 5 :
	sys.stderr.write("ERROR: missing parameter\n")
	usage()
	sys.exit(1)

pheno,threshold,perc_samples,snr_threshold,phenout = sys.argv[1:]

threshold = float(threshold)
snr_threshold = int(snr_threshold)
perc_samples = float(perc_samples)
#check files/arguments exist
if os.path.isfile(pheno) != True:
	sys.stderr.write("ERROR: file "+pheno+" not found\n")
	usage()
	sys.exit(1)

pheno = pd.read_csv(pheno, sep='\t', index_col=[0])
threshold = float(threshold)
perc_samples = float(perc_samples)

#get number of genes and samples and matrix of gene expression
samples = pheno.shape[1]
#apply filters
samples_threshold=(perc_samples/100.0)*samples
vector=SP.sum((pheno.values[:]>threshold), axis=1) #strictly greater than the gene expression threshold
booleanvector=(vector>=samples_threshold)
indexes=SP.where(booleanvector == True)[0]

print "{0} genes retained after filtering".format(len(indexes))

#filter original matrix of gene expression
filtered_matrix = pheno.iloc[indexes,:]

#select the genes based on snr threshold
index = SNR(filtered_matrix.T.values[:],snr_threshold)
final_matrix = filtered_matrix.T.iloc[:,index]
# transpose back
final_matrix = final_matrix.T[:]
#output the matrix of gene expression
final_matrix.to_csv(phenout, sep='\t', index=True)



sys.exit(0)
