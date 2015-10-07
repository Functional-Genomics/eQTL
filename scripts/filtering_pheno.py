#!/usr/bin/env python

import sys, os
import h5py
import scipy as SP
import normalise as NM


def usage():
	print """ 
This script filters the matrix of gene expression values based on a gene expression threshold (e.g. FPKMs) and minimum percentage n of samples meeting this threshold. Optionally,
it transforms filtered gene counts to gauss distribution. 

Usage:
filtering_pheno.py <pheno.hdf5> <min_expr> <min_perc_samples> <expr_transform> <pheno.filtered.hdf5> 

- min_exp = [INT | FLOAT]
- min_perc_samples = [INT | FLOAT]
- expr_transform = [gauss | log | none]"""

#check arguments
if len(sys.argv[1:]) < 5 :
	sys.stderr.write("ERROR: missing parameter\n")
	usage()
	sys.exit(1)

#arguments
pheno=sys.argv[1]
threshold=sys.argv[2]
perc_samples=sys.argv[3]
expr_transform = sys.argv[4]
phenout=sys.argv[5]  

if os.path.isfile(pheno) != True:
	sys.stderr.write("ERROR: file "+pheno+" not found\n")
	usage()
	sys.exit(1)
elif expr_transform not in [ 'gauss', 'log', 'none' ]:
        sys.stderr.write("ERROR: method "+expr_transform+" not available\n")
        usage()
        sys.exit(1)


pheno = h5py.File(pheno, 'r')
threshold = float(threshold)
perc_samples = float(perc_samples)
phenout = h5py.File(phenout, 'w')


#get number of genes and samples and matrix of gene expression
genes=pheno['phenotype/col_header/phenotype_ID'].shape[0] #catch warnings here
genesvector = pheno['phenotype/col_header/phenotype_ID'][:] #catch warnings here
samples=pheno['phenotype/row_header/sample_ID'].shape[0] #catch warnings here
samplesvector = pheno['phenotype/row_header/sample_ID'][:] #catch warnings here
matrix = pheno['phenotype/matrix'][:] #catch warnings here
chrvector = pheno['phenotype/chrom'][:] #catch warnings here
endvector = pheno['phenotype/end'][:] #catch warnings here
startvector = pheno['phenotype/start'][:] #catch warnings here

#apply filters
samples_threshold=(perc_samples*samples)/100
vector=SP.sum((matrix>=threshold), axis=0)
booleanvector=(vector>=samples_threshold)
indexes=SP.where(booleanvector == True)[0]

print "{0} genes retained after filtering".format(len(indexes))

#filter original matrix of gene expression
filtered_matrix = matrix[:,indexes]
matrix_shape=filtered_matrix.shape
filtered_genes=SP.take(genesvector,indexes, axis=0)
genes_shape=filtered_genes.shape
filtered_chrom=SP.take(chrvector, indexes, axis=0)
chrom_shape=filtered_chrom.shape
filtered_end=SP.take(endvector, indexes, axis=0)
end_shape=filtered_end.shape
filtered_start=SP.take(startvector, indexes, axis=0)
start_shape=filtered_start.shape

#apply transformation to gene expression data
if expr_transform == 'gauss':
	Y=NM.gaussianize(filtered_matrix[:])
elif expr_transform == 'log':
	Y=NM.logtransform(filtered_matrix[:])
else:
	Y=filtered_matrix[:]

#populating the output file with filtered matrix and related annotations
#filtered matrix (not-transformed)
filtered_hd_dset=phenout.create_dataset('phenotype/matrix', matrix_shape, dtype="float64")
filtered_hd_dset[...]=filtered_matrix
#tranformed values
filtered_hd_normalised=phenout.create_dataset('phenotype/Ytransformed', matrix_shape, dtype="float64")
filtered_hd_normalised[...]=Y
#samples
filtered_hd_samples=phenout.create_dataset('phenotype/row_header/sample_ID', samplesvector.shape, dtype="S100")
filtered_hd_samples[...]=samplesvector 
#genes
filtered_dset1=phenout.create_dataset('phenotype/col_header/phenotype_ID', genes_shape, dtype="S100")
filtered_dset1[...]=filtered_genes
#chrom
filtered_dset2=phenout.create_dataset('phenotype/chrom',chrom_shape, dtype="S100")
filtered_dset2[...]=filtered_chrom
#end
filtered_dset3=phenout.create_dataset('phenotype/end', end_shape, dtype="int64")
filtered_dset3[...]=filtered_end
#start
filtered_dset4=phenout.create_dataset('phenotype/start', start_shape, dtype="int64")
filtered_dset4[...]=filtered_start

phenout.close()

sys.exit(0)
