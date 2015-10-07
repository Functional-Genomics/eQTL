#!/usr/bin/env python


import sys, os
import scipy as sp
import h5py


def usage():
	print '''
This script generates an hdf5 file with a matrix of N samples x K known covariates from a csv file. Row names = sample names are required in the csv file.

Usage:

get_covariates.py <covariates.csv> <covariates.hdf5>
'''


csv,covout=sys.argv[1:]

if len(sys.argv[1:]) < 2:
	sys.stderr.write('ERROR: missing parameters\n')
	usage()
	sys.exit(1)


if os.path.isfile(csv) != True:
	sys.stderr.write('ERROR: '+csv+' not found\n')
	sys.exit(1)


#load the csv file

csvin = sp.loadtxt(csv, delimiter='\t', dtype='S100')

#open out file and write matrix of covariates and sample IDs
csvout = h5py.File(covout, 'w')
dset = csvout.create_dataset('covariates', csvin[:,1:].shape, dtype='float64')
dset[...]=csvin[:,1:].astype('float64')

dset2 = csvout.create_dataset('row_header/sample_ID', (csvin[:].shape[0],), dtype='S100')
dset2[...]=csvin[:,0][:]

csvout.close()



