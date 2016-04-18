#!/urs/bin/env python

import sys
import pandas as pd
import scipy as sp

def usage():
	print '''
This script reads a table of genotypes with variant names or positions on the first column and compare them with annotations provided in a bed4 file. Only common variants are retained and wrote to the output file.

Usage:
filter_genotype_metrix.py <genotype.tsv> <chr_pos.bed> <outfile.tsv>
'''

if len(sys.argv[1:]) < 3:
	sys.stderr.write('\n\nERROR: Missing argument\n')
	usage()
	sys.exit(1)

#reading arguments
matrixfile,bedfile,outfile = sys.argv[1:]
matrix = pd.read_csv(matrixfile,sep='\t',chunksize=10000) #read the file in chunks. Suitable for big tables
bed = pd.read_csv(bedfile,sep='\t',header=None) #read the bed file

#writing output
header = pd.read_csv(matrixfile,sep='\t',nrows=1,header=None) #read only the first line 
header = "\t".join(header.values[0].tolist())+"\n" #write header with column names: pos + samples
out = open(outfile,'w') #open outfile
out.write(header) #write header into outfile

#iterate through the file to get only positions present into the bed file
sys.stdout.write('Reading positions from the genotype matrix and comparing with the bed file...\n')
for c in matrix:
	bv = sp.in1d(c[c.columns[0]],bed[bed.columns[3]]) #take indexes of the rows to be kept based on bed file
	pos_kept = c[bv] #keep only positions in the bed file
	pos_kept.to_csv(out,mode='a',header=None,sep='\t',index=None)

sys.stdout.write('Done.')

out.close()
sys.exit(0)



	

