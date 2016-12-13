#!/usr/bin/env python


import sys, os
import h5py
import scipy as sp


def usage():
	print '''
This script sorts the samples in the pheno and covariates files based on the samples order of the geno matrix. This needs to be done before running peer/panama in the step 3 of the pipeline!

Usage:

sort_ids.py <reference.hdf5> <dataset_ref> <file.hdf5> <file.out.hdf5> <tag> 

tag: [phenotype | s_covariates | p_covariates ]
'''

if len(sys.argv[1:]) < 5:
	sys.stderr.write('ERROR: missing parameters\n')
	usage()
	sys.exit(1)

#set arguments
ref,dataset_ref,filename,filenameout,tag = sys.argv[1:]

if os.path.isfile(filename)!=True:
	sys.stderr.write('ERROR: '+filename+' not found\n')
	sys.exit(1)

if os.path.isfile(ref) != True:
	sys.sderr.write('ERROR: '+ref+' not found\n')
	sys.exit(1)

if tag not in ['phenotype','s_covariates', 'p_covariates']:
	sys.stderr.write('ERROR: unrecognised tag.\n')
	usage()
	sys.exit(1)

#open files
f= h5py.File(filename,'r')
geno_ids = h5py.File(ref, 'r')
genoids = geno_ids[dataset_ref][:]

if tag == 'phenotype':
	phenoids = f['phenotype/row_header/sample_ID'][:]
else:
	covids = f['row_header/sample_ID'][:] #for covariates

#take indices of geno ids in pheno or covariates samples
try:
	if tag == 'phenotype':
		ip = map(lambda x:(phenoids.tolist()).index(x), genoids.tolist())
	else:
		ic = map(lambda x:(covids.tolist()).index(x), genoids.tolist()) #for covariates
except:
	e = sys.exc_info()[1]
	sys.stderr.write('ERROR : {0}\n'.format(e))
	sys.exit(1)

#check all ids are included in target files
if tag == 'phenotype':
	if len(ip) != phenoids.shape[0]:
		bv = sp.in1d(phenoids,genoids)
		for x in xrange(len(phenoids[~bv])):
			sys.stderr.write('The following genotype id {0} was not found in the phenotype file\n'.format(phenoids[~bv][x]))
		sys.exit(1)
else: #for covariates
	if len(ic) != covids.shape[0]:
		bv = sp.in1d(covids,genoids)
		for x in xrange(len(covids[~bv])):
			sys.stderr.write('The following genotype id {0} was not found in the covariates file\n'.format(covids[~bv][x]))
		sys.exit(1)

#open out file
out = h5py.File(filenameout,'w')
#sort pheno or covariates ids based of geno indices and write into an out fil
if tag == 'phenotype':
	dset = out.create_dataset('phenotype/row_header/sample_ID',data=f['phenotype/row_header/sample_ID'][:][ip])
	dset = out.create_dataset('phenotype/matrix',data=f['phenotype/matrix'][:][ip,:])
	dset = out.create_dataset('phenotype/col_header/phenotype_ID',data=f['phenotype/col_header/phenotype_ID'][:])
	dset = out.create_dataset('phenotype/chrom',data=f['phenotype/chrom'][:])
	dset = out.create_dataset('phenotype/end',data=f['phenotype/end'][:])
	dset = out.create_dataset('phenotype/start',data=f['phenotype/start'][:])
else: #for all other  covariates
	dset = out.create_dataset('row_header/sample_ID',data=f['row_header/sample_ID'][:][ic])
	dset = out.create_dataset('covariates',data = f['covariates'][:][ic,:])
	dset = out.create_dataset('col_header/phenotype_ID',data=f['col_header/phenotype_ID'][:])
	#covariates['row_header/sample_ID'][:] = covariates['row_header/sample_ID'][:][ic]
	#covariates['covariates'][:] = covariates['covariates'][:][ic,:]

#close files
out.close()
sys.exit(0)

