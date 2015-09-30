#!/usr/bin/env python
import getopt, sys, os
import h5py
import pandas as ps

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


### reading (Pheno reading should be reviewed)
mapfile=np.loadtxt(mapfile, delimiter='\t', dtype='S50')

GE = pd.read_csv(phenofile, delimiter = "\t")

header = GE.colums.values.tolist()
#header = map(lambda x:x.replace('"',''), header)
rnafiles=mapfile[1:,1].tolist()
n=map(lambda x:header.index(x),rnafiles)
n.insert(0,0)

#grep only columns from the GEarray matching the RNA samples with a correspondent VCF analysis IDs.
GEsliced=GE.take(n, axis=1)

#get shape of the matrix of genes
annot=GEliced.shape[0]-1
nsamples=GEsliced.shape[1]-1
print 'matrix with {0} genes and {1} samples'.format(annot,nsamples)

#transpose the matrix
GEt = GEsliced.T
#take indexes to sort the gene expression table based on VCF samples order
#indexes_again=map(lambda x:GEarraysliced[0].tolist().index(x), samples)
#GEarraysliced[:,range(1,nsamples+1)]=GEarraysliced[:,indexes_again]

#write the new CSV file with Gene expression 
phenofileout = open(phenofileout, 'w')
Get.to_csv(phenofileout, sep = '\t', header = False)
phenofileout.close()


