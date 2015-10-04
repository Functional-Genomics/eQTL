#!/usr/bin/env python 

import sys,os 

def usage():
	print '''
This script creates a dictionary of parameters used for eQTL analysis.

Usage:
eqtl_settings.py <geno.hdf5> <pheno.filtered.hdf5> <correction_method> <correction_method.hdf5> <Kpop.hdf5> <covariates.hdf5> 
'''

def read_args(geno,pheno,correction_method,hdf5_correction,Kpop,covariates):
	CFG = {}
	#data files
	CFG['data'] = {}
	CFG['data']['geno'] = geno # chr1.hdf5
	CFG['data']['pheno']  = pheno #pheno.filtered.hdf5
	CFG['data']['kinship'] = Kpop # Kpop.hdf5
	CFG['data']['covariates'] = covariates #covariates.hdf5 (TODO: script to get known covariates)
	if correction_method == 'peer':
		CFG['data']['correction'] = hdf5_correction # peer.hdf5 which contains residuals
	elif correction_method == 'panama':
		CFG['data']['correction'] = hdf5_correction # panama.hdf5 which contains Ktot
	else:
		CFG['data']['correction'] = hdf5_correction #none.hdf5 which contains Kpop (but maybe won't be used)
	return CFG


if __name__ == '__main__':
	#the following lines are used only when the script is executed to check arguments. /
	#TODO: I would suggest to run this script before running the next eQTL scripts in the make file to check if files are ok./
	#if files are ok, then this script is imported as module in the next step.
	if len(sys.argv[1:]) < 6:
		usage()
		sys.stderr.write("ERROR: missing parameters\n")
		sys.exit(1)
	
	geno,pheno,correction_method,hdf5_correction,Kpop,covariates = sys.argv[1:]
	
	if os.path.isfile(geno) != True:
		sys.stderr.write("ERROR: "+geno+" not found\n")
		sys.exit(1)
	if os.path.isfile(pheno) != True:
                sys.stderr.write("ERROR: "+pheno+" not found\n")
                sys.exit(1)		 
        if os.path.isfile(hdf5_correction) != True:
                sys.stderr.write("ERROR: "+hdf5_correction+" not found\n")
                sys.exit(1)
        if os.path.isfile(covariates) != True:
                sys.stderr.write("ERROR: "+covariates+" not found\n")
                sys.exit(1)

	print 'Check completed.'
	sys.exit(0)


