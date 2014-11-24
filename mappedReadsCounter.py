'''
Designs oligos for a pre RNA-seq selection method
'''

### imports ###
import sys
import os
import matplotlib
import pickle
import numpy as np

		



# read in arguments

alignmentDirectoryPath = sys.argv[1] # path to a directory containing bowtie alignments and log files
backgroundFrequencyDictPath = sys.argv[2] # path to a directory containing fasta files containing sequences we don't want the oligos to hybridize to

# read in sequences
backgroundFrequencyDict = {}
for alignmentFile in os.listdir(alignmentDirectoryPath):
	if "align" in alignmentFile:
		with open(alignmentDirectoryPath + "/" + alignmentFile) as f:
			data = f.readlines()
		oligo = alignmentFile.replace("_align","")	
		#print(oligo, len(data) - 1)
		backgroundFrequencyDict[oligo] = len(data) - 1 # each oligo is expected to align once against the target genome

# read in background frequency dict
outFile = open(backgroundFrequencyDictPath, "wb")
backgroundFrequencyDict = pickle.dump(backgroundFrequencyDict, outFile)
outFile.close()

