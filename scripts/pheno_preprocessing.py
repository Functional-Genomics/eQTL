#!/usr/bin/env python
import getopt, sys, os
import h5py 
import numpy as np
import pandas as pd
import os.path
import warnings

def usage():
	print """
This script generates a matrix of gene expression with samples specified in the map file. It requires 3 mandatory arguments.
		
Usage: pheno_preprocessing.py <map_file.tsv> <pheno.tsv> <outphenoname.tsv>
		"""
#check arguments
if len(sys.argv[1:])!=3:
        usage()
        sys.exit(1)

#arguments
mapfile=sys.argv[1]
phenofile=sys.argv[2]
phenofileout=sys.argv[3]

if os.path.isfile(mapfile)!=True:
        sys.stderr.write("ERROR: mapfile "+mapfile+" not found\n")
        sys.exit(1)

if os.path.isfile(phenofile)!=True:
        sys.stderr.write("ERROR: phenofile "+phenofile+" not found\n")
        sys.exit(1)


##############################################
### reading (Pheno reading should be reviewed)
try:
        with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                mapfile=np.loadtxt(mapfile, delimiter='\t', dtype='S50')
except Exception as e:
        sys.stderr.write("ERROR: unable to load "+mapfile+" - "+str(e)+"\n")
        sys.exit(1)

try:
        GE = pd.read_csv(phenofile, delimiter = "\t")
except Exception as e:
        sys.stderr.write("ERROR: unable to load "+phenofile+" - "+str(e)+"\n")
        sys.exit(1)

header = GE.columns.values.tolist()
#header = map(lambda x:x.replace('"',''), header)
rnafiles=mapfile[1:,1].tolist()
#sanity check
if len(rnafiles) > len(header):
	sys.stderr.write("ERROR: inconsistency between number of samples in mapfile and gene expression matrix\n")
	sys.exit(1)

#take indexes of mapfile RNA samples within gene expression matrix header 
n=map(lambda x:header.index(x),rnafiles)
####
#insert index 0 in list to take also first column name
n.insert(0,0)
#grep only columns from the GEarray matching the RNA samples with a correspondent VCF analysis IDs.
GEsliced=GE.take(n, axis=1)
#susbstitute RNA samples ID with DNA samples ID
GEsliced.columns=mapfile[:,0]
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
GEt.to_csv(phenofileout, sep = '\t', header = False)
phenofileout.close()
sys.exit(0)
