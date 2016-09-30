#!/usr/bin/env python

import fileinput
import sys, os

def usage():
	print """
This script filters by VAF and LOWSUPPORT and adds one new column to the original vcf file with the VAF value and print to stdout

Usage:
filter_somatic_merged_vcf.py <file.vcf.gz> <VAF_threshold> <chr> <use_pass> > filtered.vcf
use_pass= [y | n]
"""


if __name__ == "__main__":
	if len(sys.argv[1:]) < 4:
		sys.stderr.write('\nERROR: missing argument')
		usage()
		sys.exit(1)
	if  os.path.isfile(sys.argv[1]) != True:
		sys.stderr.write('ERROR: file '+sys.argv[1]+' not found\n')
		sys.exit(1)
	
	try:
		 threshold=float(sys.argv[2])
	except:
		sys.stderr.write('ERROR: please provide an integer/float value as VAF threshold\n')
		sys.exit(1)
	
	chromosome=str(sys.argv[3])
	try:
		use_pass=str(sys.argv[4])
	except:
		sys.stderr.write('ERROR: please provide y or n for use_pass parameter\n')
		sys.exit(1)

	for line in fileinput.input(sys.argv[1],openhook=fileinput.hook_compressed):
		if line.startswith('#'):
			print line.strip()
		else:
			lines=line.split('\t')
			chr=lines[0]
			if str(chr) != chromosome:
				continue
			else:
				filter=lines[-2]
				if use_pass == 'y':
					if filter == 'LOWSUPPORT':
						continue
				else:
					pass
				info=lines[-1]
				info=info.split(';')
				for v in info:
					if 'VAF' in v:
						AF = float(v.split('=')[-1])
						if AF>= threshold:
							new_list=line.strip()
							sys.stdout.write(new_list+'\n')
						else:
							continue
			
	sys.exit(0)	
