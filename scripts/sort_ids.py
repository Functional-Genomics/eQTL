#!/usr/bin/env python


import sys, os
import h5py
import scipy as sp


def usage():
	print '''
This script sorts the samples in the pheno and covariates files based on the samples order of the geno matrix. This needs to be done before running peer/panama in the step 3 of the pipeline!

Usage:

sort_ids.py <pheno.filtered.hdf5> <covariates.hdf5> <samples.hdf5> 

'''

if len(sys.argv[1:]) < 3:
	sys.stderr.write('ERROR: missing parameters\n')
	usage()
	sys.exit(1)


phenofile,covfile,samples = sys.argv[1:]

if os.path.isfile(phenofile) != True:
	sys.stderr.write('ERROR: '+phenofile+' not found\n')
	sys.exit(1)

if os.path.isfile(covfile) != True:
        sys.stderr.write('ERROR: '+covfile+' not found\n')
        sys.exit(1)

if os.path.isfile(samples) != True:
	sys.sderr.write('ERROR: '+samples+' not found\n')
	sys.exit(1)

#open files
pheno = h5py.File(phenofile)
covariates = h5py.File(covfile)
geno_ids = h5py.File(samples, 'r')

phenoids = pheno['phenotype/row_header/sample_ID'][:]
covids = covariates['row_header/sample_ID'][:]
genoids = geno_ids['sample_ID'][:]

#take indices of geno ids in pheno and covariates samples
try:
	ip = map(lambda x:(phenoids.tolist()).index(x), genoids.tolist())
	ic = map(lambda x:(covids.tolist()).index(x), genoids.tolist())
except:
	e = sys.exc_info()[1]
	sys.stderr.write('ERROR : {0}\n'.format(e))
	sys.exit(1)

if len(ip) != phenoids.shape[0]:
	bv = sp.in1d(phenoids,genoids)
	for x in xrange(len(phenoids[~bv])):
		sys.stderr.write('The following genotype id {0} was not found in the phenotype file\n'.format(phenoids[~bv][x]))
	sys.exit(1)

if len(ic) != covids.shape[0]:
	bv = sp.in1d(covids,genoids)
	for x in xrange(len(phenoids[~bv])):
		sys.stderr.write('The following genotype id {0} was not found in the phenotype file\n'.format(phenoids[~bv][x]))
	sys.exit(1)

#sort pheno and covariates ids based of geno indices
pheno['phenotype/row_header/sample_ID'][:]=pheno['phenotype/row_header/sample_ID'][:][ip]
pheno['phenotype/matrix'][:] = pheno['phenotype/matrix'][:][ip,:]
#pheno['phenotype/Ytransformed'][:] = pheno['phenotype/Ytransformed'][:][ip,:]
covariates['row_header/sample_ID'][:] = covariates['row_header/sample_ID'][:][ic]
covariates['covariates'][:] = covariates['covariates'][:][ic,:]

#close files
pheno.close()
covariates.close()

sys.exit(0)

