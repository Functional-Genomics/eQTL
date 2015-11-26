#!/usr/bin/env python

import fileinput
#import sys, os

def usage():
	print ''' 
This script substitutes the uncalled variants with "0" or "0/0" and low quality variants with "." in the merged vcf file.

Usage:

substitute_vcf.py <chr1.merged.vcf> '''

def lowqual2missing(genotype):
	if genotype.strip()=='0':
		genotype='.'
	else:
		pass
	return genotype

if __name__ == "__main__":

#	if len(sys.argv[1:])<1:
#		usage()
#		sys.stderr.write("\nERROR: missing file\n")
#		sys.exit(1)
#
#	vcf=sys.argv[1] 

#	if os.path.isfile(vcf) != True:
#		sys.stderr.write("\nERROR: file "+vcf+" not found\n")
#		sys.exit(1)


	#with open(sys.stdin,'r') as f:
	#	for line in f:
	for line in fileinput.input():
		if line.startswith('#'):
			print line.strip()
		else:
			line=line.split('\t')
			o=map(lambda x:x.strip().replace('.','-') or x.strip().replace('./.','-/-'),line[9:])
			o=map(lambda x:lowqual2missing(x),o)
			o=map(lambda x:x.strip().replace('-','0'),o) #replace uncalled with reference allele
			o="\t".join(o)
			a="\t".join(line[:9])
			print a+'\t'+o


	sys.exit(0)	
