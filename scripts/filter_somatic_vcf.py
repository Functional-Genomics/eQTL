#!/usr/bin/env python

import fileinput
import sys, os

def usage():
	print """
This script filters by VAF and LOWSUPPORT and adds one new column to the original vcf file with the VAF value. 

Usage:
filter_somatic_merged_vcf.py <file.vcf> <VAF_threshold <chr> > filtered.vcf
"""


if __name__ == "__main__":
	if len(sys.argv[1:]) < 3:
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
				if filter == 'LOWSUPPORT':
					continue
				else:
					info=lines[-1]
					info=info.split(';')
					for v in info:
						if 'VAF' in v:
							AF = float(v.split('=')[-1])
							if AF>= threshold:
								new_list=line.strip()+'\t'+str(AF)
								print new_list
							else:
								continue
			
	sys.exit(0)	
