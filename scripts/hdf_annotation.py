#!/usr/bin/env python 

import sys, os
import h5py
import scipy as sp

#check number of arguments
if len(sys.argv[1:]) != 2 :
	sys.stderr.write("ERROR: full set of parameters not provided \n")
	usage()
	sys.exit(1)


#arguments
annotation=sys.argv[1]
phenofile=sys.argv[2]

#check files exist
if os.path.isfile(annotation) != True
	sys.stderr.write("ERROR: file "+annotation+" not found\n")
	sys.exit(1)
elif os.path.isfile(phenofile) != True:
	sys.stderr.write("ERROR: file "+phenofile+" not found\n")
	sys.exit(1)

#open files
annotation = sp.loadtxt(annotation, delimiter='\t', dtype='S100')
phenofile=h5py.File(phenofile) #appending mode

#take annotations for genes in gene matrix
booleanvector=sp.in1d(annotation[:,3][:], phenofile['phenotype/col_header/phenotype_ID'][:].tolist()) #catch warnings if hdf5 dataset does not exist
#chrom
dataset=phenofile.create_dataset('phenotype/chrom', (sum(booleanvector),), dtype='S10') #same here
dataset[...]=annotation[sp.where(booleanvector)[0],0]
#start
dataset=phenofile.create_dataset('phenotype/start', (sum(booleanvector),), dtype='int64') #same here
dataset[...]=annotation[sp.where(booleanvector)[0],1]
#end
dataset=phenofile.create_dataset('phenotype/end', (sum(booleanvector),), dtype='int64') #same here
dataset[...]=annotation[sp.where(booleanvector)[0],2]

#check keys (maybe useless)
for i in phenofile.keys():
	for x in phenofile[i].keys():
		print x

phenofile.close()


