'''
Designs oligos for a pre RNA-seq selection method
'''

### imports ###
import sys
import os
import matplotlib
import pickle
import numpy as np
matplotlib.use('Agg')
from matplotlib import pyplot as plt


def readFastaFile(fastaFilePath):
	'''
	Given a path to a multiline fasta file, reads the file, returning two lists - one containing the sequences, the other containing the headers
	inputs: path to a fasta file
	outputs: a list of the sequences, a list of the sequence headers
	'''
	sequences = []
	headers = []
	with open(fastaFilePath) as f:
		data = f.readlines()
	sequence = ""
	for line in data:
		if ">" in line:
			header = line.replace(">", "").strip()
			headers.append(header)
			if not sequence == "":
				sequences.append(sequence.upper())
				sequence = ""
		else:
			sequence += line.strip()
	sequences.append(sequence.upper())

	return sequences, headers

def designOligos(targetSequences, backgroundFrequencyDict, targetLength, maxGapLength, threshold, outputPath):
	'''
	Given the frequency of oligos in the backgroundFrequencyDict, construct oligos of length targetLength that tile the input sequences with a maximum distance of maxGapLength between them
	inputs: dictionary that gives the frequency of all background oligos, target sequences we want to design oligos to hybridize to, target length of the oligos, max gap length between the oligos
	outputs: writes the designed oligos to a Fasta file

	'''
	designedOligos = []
	designedPositions= []
	backgroundCounts = list(backgroundFrequencyDict.values())
	designedOligoFrequencies = [] # how many times do the designed oligos appear in the background

	seenOligos = set()

	for i in range(len(targetSequences)):
		currentSequence = targetSequences[i]
		pos = 0
		maxIndex = len(currentSequence) - targetLength - 1

		# iterate until sequence is fully tiled; maximizes density while avoiding overlapping oligos
		while pos <= maxIndex:
			# generate all oligo options for this window
			j = 0
			oligoTuples = []
			while j < maxGapLength and pos + j <= maxIndex:
				oligo = currentSequence[pos + j  : pos + j + targetLength] 
				if not oligo in seenOligos: # only use each oligo once
					if oligo in backgroundFrequencyDict:
						freq = backgroundFrequencyDict[oligo]
					else:
						freq = 0
					if freq < threshold:
						oligoPos = pos + j + targetLength
						oligoTuples.append((freq, oligoPos, oligo))
				j += 1
			# pick the oligo with the lowest frequency
			if oligoTuples:
				oligoTuples.sort()
				selectedOligoTuple = oligoTuples[0]
				pos = selectedOligoTuple[1]	
				designedPositions.append(pos)
				designedOligos.append(selectedOligoTuple[2])
				designedOligoFrequencies.append(selectedOligoTuple[0])
				seenOligos.add(selectedOligoTuple[2])
			else:
				# no oligos satisfy design requirements
				pos = pos + maxGapLength + targetLength
	# write fasta file
	outFile = open(outputPath + "/designedOligos.fa", "w")
	complementFile = open(outputPath + "/coveredRegions.fa", "w")
	for i in range(len(designedOligos)):
		outFile.write(">oligo_" + str(i) + "\t" + str(designedOligoFrequencies[i]) + "\t" + str(designedPositions[i])+ "\n")
		outFile.write(revComp(designedOligos[i]) + "\n")
		complementFile.write(">oligo_" + str(i) + "\t" + str(designedOligoFrequencies[i]) + "\t" + str(designedPositions[i])+ "\n")
		complementFile.write(designedOligos[i] + "\n")
	outFile.close()
	complementFile.close()
	return designedOligoFrequencies, designedPositions

def revComp(seq):
	'''
	reverse complements a sequence
	inputs: a sequence
	outputs: a reverse compelemented sequence
	'''
	rcDict = {"A":"T", "T":"A", "U":"A", "C":"G", "G":"C"}
	toReturn = ""
	for char in seq:
		toReturn = rcDict[char] + toReturn
	return toReturn
		



if __name__ == "__main__":
	targetDirectoryPath = sys.argv[1] # path to a directory containing fasta files giving the sequences we want the oligos to hybridize to
	backgroundFrequencyDictPath = sys.argv[2] # path to a directory containing fasta files containing sequences we don't want the oligos to hybridize to
	targetOligoLength = int(sys.argv[3]) # desired length of oligos
	maxGapLength = int(sys.argv[4]) # desired length of oligos
	threshold = int(sys.argv[5])
	outputPath = sys.argv[6] # path to write output files
	
	# create output directory if it doesn't exist
	if not os.path.isdir(outputPath):
		os.makedirs(outputPath)

	# intialize lists
	allTargetSequences = []
	allTargetHeaders = []
	allBackgroundSequences = []
	allBackgroundHeaders = []

	# read in sequences
	print("reading target files")
	for targetFile in os.listdir(targetDirectoryPath):
		print(targetFile)	
		targetSequences, targetHeaders = readFastaFile(targetDirectoryPath + "/" + targetFile)
		allTargetSequences += targetSequences
		allTargetHeaders += targetHeaders
	
	# read in background frequency dict
	backgroundFrequencyDict = pickle.load(open(backgroundFrequencyDictPath, "rb"))

	print("designing oligos")	
	designedOligoFrequencies, designedPositions = designOligos(allTargetSequences, backgroundFrequencyDict, targetOligoLength, maxGapLength, threshold, outputPath)
		
	print("making plots")
	# plot background frequencies
	plt.hist(list(backgroundFrequencyDict.values()))
	plt.xlabel("Oligo Frequency")
	plt.ylabel("Frequency")
	plt.savefig(outputPath + "/backgroundFrequency.png")
	plt.close()

	# plot oligo frequency distribution vs background
	plt.hist(designedOligoFrequencies)
	plt.title("Hybridization to Background (Mean="+str(np.mean(designedOligoFrequencies))+ ", Std="+str(np.std(designedOligoFrequencies)) + ")")
	plt.xlabel("Oligo Frequency")
	plt.ylabel("Frequency")
	plt.savefig(outputPath + "/designedOligoFrequency.png")
	plt.close()

	# plot coverage
	# calculate differences
	diff = [designedPositions[n] - designedPositions[n-1] for n in range(1, len(designedPositions))]
	plt.hist(diff)
	plt.title("Coverage (Mean="+str(np.mean(diff))+ ", Std="+str(np.std(diff)) + ")")
	plt.xlabel("Distance between adjacent oligo")
	plt.ylabel("Frequency")
	plt.savefig(outputPath + "/coverage.png")
	plt.close()
	


