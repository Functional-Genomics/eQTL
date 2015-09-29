#!/usr/bin/env python
import getopt, sys, os
import h5py
import csv
import numpy as np

def usage():
	print """
		 This script generates a csv matrix of gene expression with samples specified in the map file. It requires 3 mandatory arguments.
		
		Usage:
		pheno_preprocessing.py <map_file.csv> <pheno.tsv> <outphenoname.matched.csv>
		"""

if len(sys.argv[1:])!=3:
        usage()
        sys.exit()

#arguments
mapfile=sys.argv[1]
phenofile=sys.argv[2]
phenofileout=sys.argv[3]


### reading
mapfile=np.loadtxt(mapfile, delimiter='\t', dtype='S50')

GE=open(file,'r').readlines()
GE=[i.strip().split('\t') for i in GE]
GE[0]=map(lambda x:x.replace('"',''), GE[0])

rnafiles=mapfile[1:,1].tolist()
n=map(lambda x:GE[0].index(x),rnafiles)
n.insert(0,0)

#grep only columns from the GEarray matching the RNA samples with a correspondent VCF analysis IDs.
GEarray=np.array(GE)
GEarraysliced=np.take(GEarray, n, axis=1)
annot=GEarraysliced.shape[0]
nsamples=GEarraysliced.shape[1]-1

print 'matrix with {0} genes and {1} samples'.format(annot,nsamples)

#take indexes to sort the gene expression table based on VCF samples order
#indexes_again=map(lambda x:GEarraysliced[0].tolist().index(x), samples)
#GEarraysliced[:,range(1,nsamples+1)]=GEarraysliced[:,indexes_again]
GEt=GEarraysliced.transpose()
GEt[0][0]='Samples'
#write the new CSV file with Gene expression 
with open(phenofileout, 'w') as csvfile:
    GEwriter = csv.writer(csvfile, delimiter='\t', quotechar='|',quoting=csv.QUOTE_MINIMAL)
    for row in GEt:
        GEwriter.writerow(row)

phenofileout.close()


