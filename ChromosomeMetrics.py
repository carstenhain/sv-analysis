import numpy as np

#Returns all canonical human chromosomes as array
def GetChromosomes ():	
	return ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

#Returns all canonical human chromosomes minus chrX and chrY as array
def GetChromosomesWOXY ():	
	return ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]

#Returns all canonical human chromosome with their length as dictionary
def GetChromosomeLengths ():	
	return { 
		"chr1" : 248956422,
		"chr2" : 242193529,
		"chr3" : 198295559,
		"chr4" : 190214555,
		"chr5" : 181538259,
		"chr6" : 170805979,
		"chr7" : 159345973,
		"chr8" : 145138636,
		"chr9" : 138394717,
		"chr10" : 133797422,
		"chr11" : 135086622,
		"chr12" : 133275309,
		"chr13" : 114364328,
		"chr14" : 107043718,
		"chr15" : 101991189,
		"chr16" : 90338345,
		"chr17" : 83257441,
		"chr18" : 80373285,
		"chr19" : 58617616,
		"chr20" : 64444167,
		"chr21" : 46709983,
		"chr22" : 50818468,
		"chrX" : 156040895,
		"chrY" : 57227415
	}

def GetGenomeLength ():
	return 3088269832

#Returns absolute position for a chromosome + position combination as int
def GetAbsolutePosition (chrom, pos):	
	abs_pos = np.int64(pos)	
	for chromosome in GetChromosomes():		
		if chromosome == chrom:			
			break			
		abs_pos += np.int64(GetChromosomeLengths()[chromosome])		
	return abs_pos
