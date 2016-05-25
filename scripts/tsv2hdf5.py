#!/usr/bin/env python

import sys,os
import h5py
import scipy as sp
import pandas as pd

def usage():
	print '''

Usage:
tsv2hdf5.py <file.tsv> <dataset_name> <row_names> <col_names> <boolean> <outfile.hdf5>
'''

if len(sys.argv[1:])<4:
	usage()
	sys.stderr.write('\nERROR: missing argument\n')
	sys.exit(1)


infile, dataset_name, row_names,col_names, transpose, outfile = sys.argv[1:]

if os.path.isfile(infile) != True:
	sys.stderr.write('\nERROR: file '+file2+' not found\n')
	sys.exit(1)
if transpose != 'y' and transpose != 'n':
	sys.stderr.write('\nERROR: please enter valid boolean values: y/n\n')
	sys.exit(1)

sys.stderr.write('Reading the tsv file...Done\n')
infile = pd.read_csv(infile,sep="\t",index_col=[0])

if transpose == 'y':
	matrix = infile.values.T
	columns = infile.index.values.astype(str)
	rows = infile.columns.values.astype(str)
else:
	matrix = infile.values
	columns = infile.columns.values.astype(str)
	rows = infile.index.values.astype(str)

sys.stderr.write('Writing the hdf5 file...\n')
outfile = h5py.File(outfile,'w')
outfile.create_dataset(dataset_name,data=matrix)
outfile.create_dataset(row_names,data=rows)
outfile.create_dataset(col_names,data=columns)
outfile.close()

sys.stderr.write('Done\n')
sys.exit(0)

