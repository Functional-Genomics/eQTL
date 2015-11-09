#!/usr/bin/env python 

import sys, os
import h5py
import scipy as sp

def usage():
	print ''' This script loads annotations from tsv file into hdf5 matrix with gene expression values.

Usage:
hdf_annotation.py <annotation.tsv> <pheno.hdf5> '''


#check number of arguments
if len(sys.argv[1:]) != 2 :
	usage()
	sys.stderr.write("ERROR: full set of parameters not provided \n")
	sys.exit(1)


#arguments
annotation=sys.argv[1]
phenofile=sys.argv[2]

#check files exist
if os.path.isfile(annotation) != True:
	sys.stderr.write("ERROR: file "+annotation+" not found\n")
	sys.exit(1)
elif os.path.isfile(phenofile) != True:
	sys.stderr.write("ERROR: file "+phenofile+" not found\n")
	sys.exit(1)

#open files
annotation = sp.loadtxt(annotation, delimiter='\t', dtype='S100')
pheno=h5py.File(phenofile) #appending mode

#number of genes in annotation

print 'number of genes in annotation file {0}'.format(len(annotation))

#strip unwanted char from genes in annotation array
annotation[:,3] = sp.char.replace(annotation[:,3], '"','')

#take annotations for genes in gene matrix
booleanvector=sp.in1d(annotation[:,3][:], pheno['phenotype/col_header/phenotype_ID'][:].tolist()) #catch warnings if hdf5 dataset does not exist
#get shape
shape=annotation[sp.where(booleanvector)[0]][:,0].shape
#get dataframe from annotation

annotation=annotation[sp.where(booleanvector)[0]]

print 'number of genes kept from annotation {0}'.format(len(annotation))

annotation_list=annotation[:,3].tolist()
gene_list=pheno['phenotype/col_header/phenotype_ID'][:].tolist()
idx=[]
for gene in gene_list:
	i=annotation_list.index(gene)
	idx.append(i)

annotation=annotation[idx,:]

#chrom
dataset=pheno.create_dataset('phenotype/chrom', shape, dtype='S10') #same here
dataset[...]=annotation[:,0][:]
#start
start=(annotation[:,1]).astype('int64')
dataset=pheno.create_dataset('phenotype/start', shape, dtype='int64') #same here
dataset[...]=start[:]
#end
end=(annotation[:,2]).astype('int64')
dataset=pheno.create_dataset('phenotype/end', shape, dtype='int64') #same here
dataset[...]=end[:]

#check annotation in gene expression file

print 'Chr\tStart\tEnd\tGene'

for x in xrange(10):
	o='\t'.join(annotation[x])
	print o

#check keys (maybe useless)
print 'check if {0} has 6 datasets: matrix, col_header, row_header, chrom, start, end'.format(phenofile)
exp_keys=['matrix', 'col_header', 'row_header', 'chrom', 'start', 'end']
keys=[]
for i in pheno.keys():
	for x in pheno[i].keys():
		keys.append(x)

try:
	n = map(lambda x:keys.index(x), exp_keys)
	pheno.close()
except:
	pheno.close()
	sys.stderr.write('ERROR: dataset missing in the '+phenofile+' file\n')
	sys.exit(1)


