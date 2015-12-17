#!/usr/bin/env python

import sys,os
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def usage():
	print ''' 

This script generates a barplot of vcf filtering statistics with:
1) Mean + Stdv of Number of Variants before filtering
2) Mean + Stdv of Number of Variants after first filtering
3) Number of Variants filtered after setting missingness 
4) Number of Variants filtered after setting mac 
5) Number of Variants filtered applying both Missingness and Mac on the same file 

Usage:

get_barplot.py <vcf_snps_0.tsv> <vcf_snps_1.tsv> <vcf_snps_2.tsv> <outfilename>

TODO: Generate aggregated counts for each study!'''

def openfile(file):
	a = np.loadtxt(file,dtype='S1000') #the tsv file contains tab only in the header!
	return a




if __name__ == '__main__':
	
	if len(sys.argv[1:]) < 4:
		usage()
		sys.stderr.write('\nERROR: missing parameter\n')
		sys.exit(1)

	f1,f2,f3,outfile = sys.argv[1:]

	#if f1,f2,f3 does not exist then exit
	if os.path.isfile(f1) != True:
		sys.stderr.write('\nERROR: file '+f1+' not found\n')
		sys.exit(1)
        if os.path.isfile(f2) != True:
                sys.stderr.write('\nERROR: file '+f2+' not found\n')
                sys.exit(1)
        if os.path.isfile(f3) != True:
                sys.stderr.write('\nERROR: file '+f3+' not found\n')
                sys.exit(1)

	f1=openfile(f1)
	f2=openfile(f2)
	f3=openfile(f3)

	#generate dictionary with values for each chromosome.TODO: clean the code here
	dic={}
	for i in f1[1:]:
		key=int(i[0])
		dic[key]=[np.array(i[1:],dtype='float')] #create key and append first list of values

	for i in f2[1:]:
		key=int(i[0])
		if key in dic:
			dic[key].append(np.array(i[1:],dtype='float'))
		else:
			dic[key]=[np.array(0)]

	for i in f3.T[1:,:]:
		key=int(i[0])
		if key in dic:
			dic[key].extend((np.array(i[1:],dtype='float')).tolist())
		else:
			dic[key]=[np.array(0)]

	#load data in dframe		
	df = pd.DataFrame(dic)
	#get list of values for each filtering
	NOfilter = df.ix[0,0:] 
	PASSfilter = df.ix[1,0:]
	MISSfilter = (df.ix[2,0:]).tolist()
	MACfilter = (df.ix[3,0:]).tolist()
	MISSANDMACfilter = (df.ix[4,0:]).tolist()

	C= map(lambda x:x[0], NOfilter.iteritems())
	M0 = map(lambda x:x[1].mean(), NOfilter.iteritems())
	S0 = map(lambda x:x[1].std(), NOfilter.iteritems())
	M1 = map(lambda x:x[1].mean(), PASSfilter.iteritems())
	S1 = map(lambda x:x[1].std(), PASSfilter.iteritems())

	#generate plot. TODO: move this within funciotn
	n_groups = len(dic)
	plt.figure(figsize=(10,10))
	#TODO:clean up the code here a bit
	plt.subplot(111).spines['top'].set_visible(False)
	plt.subplot(111).spines['right'].set_visible(False)
	plt.tick_params(axis='x',which='both',bottom='off',top='off', labelsize=14)
	plt.tick_params(axis='y',which='both',left='on',right='off', labelsize=14)
	#font = {'family' : 'normal',
        #'weight' : 'bold',
        #'size'   : '15'}

	#rc('font', **font)  # pass in the font dict as kwargs
	#index = np.arange(n_groups)
	index = np.arange(0, n_groups * 2, 2)
	bar_width = 0.2

	opacity = 0.4
	error_config = {'ecolor': '0.3'}

	rects1 = plt.bar(index, tuple(M0), bar_width,
			 alpha=opacity,
			 color='b',
			 yerr=tuple(S0),
			 error_kw=error_config,
			 label='No filter')

	rects2 = plt.bar(index + bar_width, tuple(M1), bar_width,
			 alpha=opacity,
			 color='r',
			 yerr=tuple(S1),
			 error_kw=error_config,
			 label='Filter PASS variants')

	rects3 = plt.bar(index + (2*bar_width), tuple(MISSfilter), bar_width,
			 alpha=opacity,
			 color='g',
			 error_kw=error_config,
			 label='Filter Missing Genotypes')

	rects4 = plt.bar(index + (3*bar_width), tuple(MACfilter), bar_width,
			 alpha=opacity,
			 color='k',
			 error_kw=error_config,
			 label='Filter MAC')

	rects4 = plt.bar(index + (4*bar_width), tuple(MISSANDMACfilter), bar_width,
			 alpha=opacity,
			 color='m',
			 error_kw=error_config,
			 label='Filter Missing Genotype + MAC')

	#plt.grid(True)i
	plt.legend(loc='upper center', bbox_to_anchor=(0.5,1.1), ncol=3, fancybox=True, shadow=False, frameon=False, fontsize=12)

	plt.xlabel('Chr', fontsize=12)
	plt.ylabel('Variants', fontsize=12)
	plt.tick_params(direction='out')
	plt.xticks(index + (2.5*bar_width), tuple(C), fontsize=12)
	sizey=plt.ylim()
	plt.yticks(fontsize=12)
	#plt.tight_layout()
	plt.savefig(outfile, format='png',dpi=350)

	sys.exit(0)

