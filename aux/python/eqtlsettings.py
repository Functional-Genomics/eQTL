import sys,os 

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
	return CFG,correction_method
