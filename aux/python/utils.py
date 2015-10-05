import sys
import scipy as SP
import h5py
import pdb
import copy

def smartAppend(table,name,value):
	"""
	helper function for apppending in a dictionary	
	"""	
	if name not in table.keys():
		table[name] = []
	table[name].append(value)

def smartConcatenate(table,name,value):
	"""
	helper function for apppending in a dictionary	
	"""
	if name not in table.keys():
		table[name] = SP.zeros((value.shape[0],0))
	if len(value.shape)==1:
		value = value[:,SP.newaxis]
	table[name] = SP.concatenate((table[name],value),1)

def dumpDictHdf5(RV,o):
	""" Dump a dictionary where each page is a list or an array """
	for key in RV.keys():
		o.create_dataset(name=key,data=SP.array(RV[key]),chunks=True,compression='gzip')

def smartDumpDictHdf5(RV,o):
	""" Dump a dictionary where each page is a list or an array or still a dictionary (in this case, it iterates)"""
	for key in RV.keys():
		if type(RV[key])==dict:
			g = o.create_group(key)
			smartDumpDictHdf5(RV[key],g)
		else:
			o.create_dataset(name=key,data=SP.array(RV[key]),chunks=True,compression='gzip')

def getLambda(pv):
	"""
	return lambda genomic control given the pvs
	"""
	rv = SP.array([SP.median(SP.log10(pv),1)/SP.log10(0.5)])
	return rv

def getRelativePosition(pos,strand,start,end):
	""" return the relative position respect to the TSS """
	if strand=='1':
		rv=float(pos-start)
	if strand=='-1':
		rv=float(end-pos)
	else:
		rv = 1.009
	return rv

def matchIDs(ID1,ID2):
	""" match ID1 and ID2 """
	idx1 = []
	idx2 = []
	for p in range(ID1.shape[0]):
		is_in = ID1[p] in ID2
		idx1.append(is_in)
		if is_in:	idx2.append(SP.where(ID2==ID1[p])[0][0])
	idx1 = SP.array(idx1)
	idx2 = SP.array(idx2)
	return idx1, idx2

def wait(sec,verbose=True):
	""" wait sec seconds """
	import time as TIME
	if verbose:
		print "wait %s s"%sec
	start = TIME.time()
	while 1:
		if TIME.time()-start>sec:	break
	pass

def pearsCorrRavel(Y1,Y2):
	""" calculated the prearson correlation between vec(Y1) and vec(Y2) """

	y1 = Y1.ravel()
	y2 = Y2.ravel()
	rv = SP.corrcoef(y1,y2)[0,1]

	return rv


def pearsCorrMean(Y1,Y2):
	""" calculated the avg prearson correlation between columns of Y1 and Y2 """

	rv = 0
	for ic in range(Y.shape[1]):
		rv += SP.corrcoef(Y1[:,ic],Y2[:,ic])[0,1]
	rv/=float(Y.shape[1])

	return rv

def getCumPos(chrom,pos):
	"""
	getCumulativePosition
	"""
	n_chroms = int(chrom.max())
	x = 0
	posCum = copy.copy(pos)
	for chrom_i in range(1,n_chroms+1):
		I = chrom==chrom_i
		posCum[I]+=x
		x=posCum[I].max()
	return posCum

def getChromBounds(chrom,posCum):
	"""
	getChromBounds
	"""
	n_chroms = int(chrom.max())
	chrom_bounds = []
	for chrom_i in range(2,n_chroms+1):
		I1 = chrom==chrom_i
		I0 = chrom==chrom_i-1
		_cb = 0.5*(posCum[I0].max()+posCum[I1].min())
		chrom_bounds.append(_cb)
	return chrom_bounds


