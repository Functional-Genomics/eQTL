#!/usr/bin/env python

import sys,os
import h5py
import scipy as sp


def usage():
	print '''
This script checks if number of DNA samples in the mapfile is equal to the number of samples in the samples.hdf5 file.

Usage:

check_consistency.py <map_file.csv> <samples.hdf5>

'''


if len(sys.argv[1:]) < 2:
	sys.stderr.write('ERROR: missing paramters\n')
	usage()
	sys.exit(1)

mapfile,samples=sys.argv[1:]

if os.path.isfile(mapfile) !=True:
	sys.stderr.write('ERROR: '+mapfile+' not found\n')
	sys.exit(1)

if os.path.isfile(samples) != True:
	sys.stderr.write('ERROR: '+samples+' not found\n')
	sys.exit(1)


mapfile=sp.loadtxt(mapfile, delimiter='\t', dtype='S50')
samples=h5py.File(samples,'r')

DNA_id=mapfile[1:,0][:]
sample_id=samples['sample_ID'][:]

try:
	#check if ALL sample IDs in mapfile are in hdf5 files from chromosomes
	n=map(lambda x:(sample_id.tolist()).index(x), DNA_id.tolist())
	sys.exit(0)
except:
	#if not exit
	sys.stderr.write('ERROR: inconsistency between samples.hdf5 and mapfile\n')
	sys.exit(1) 





