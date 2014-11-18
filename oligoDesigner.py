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

def countOligos(sequences, targetLength):
	'''
	Given a list of sequences, creates a dictionary that gives the number of occurences each oligomer of the given targetLength appears
	inputs: list of sequences, targetLength of mers to count
	outputs: a dictionary of the following format - {oligomer: frequency}
	'''
	frequencyDict = {}
	for i in range(len(sequences)):
		seq = sequences[i]
		seqtargetLength = len(seq)
		if seqtargetLength >= targetLength:
			for j in range(seqtargetLength - targetLength):
				print(i, j, seqtargetLength)
				oligo = seq[j: j + targetLength]
				if oligo in frequencyDict:
					frequencyDict[oligo] += 1
				else: 
					frequencyDict[oligo] = 1

	return frequencyDict

def designOligos(targetSequences, backgroundFrequencyDict, targetLength, maxGapLength, outputPath):
	'''
	Given the frequency of oligos in the backgroundFrequencyDict, construct oligos of length targetLength that tile the input sequences with a maximum distance of maxGapLength between them
	inputs: dictionary that gives the frequency of all background oligos, target sequences we want to design oligos to hybridize to, target length of the oligos, max gap length between the oligos
	outputs: writes the designed oligos to a Fasta file

	'''
	designedOligos = []
	backgroundCounts = list(backgroundFrequencyDict.values())
	designedOligoFrequencies = [] # how many times do the designed oligos appear in the background
# not sure if a threshold is neccessarily, just use the oligo with the lowest frequency in the background
#	mean = np.mean(backgroundCounts)
#	std =  np.std(backgroundCounts)
#	threshold = mean - 3 * std # values lower than 3 standard devations of the mean should constitute just 0.0015% of all the total values
#	print("Using threshold: " + str(threshold))
#	if threshold < 0:
#		threshold = 0
	seenOligos = set()

	for i in range(len(targetSequences)):
		currentSequence = targetSequence[i]
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
					oligoPos = pos + j + targetLength
					oligoTuples.append((freq, oligoPos, oligo))
				j += 1
			# pick the oligo with the lowest frequency
			oligoTuples.sort()
			selectedOligoTuple = oligoTuples[0]
			pos = selectedOligoTuple[1]	
			designedOligos.append(selectedOligoTuple[2])
			designedOligoFrequencies.append(selectedOligoTuple[0])
			seenOligos.add(selectedOligoTuple[2])
	# write fasta file
	outFile = open(outputPath + "/designedOligos.fa", "w")
	for i in range(len(designedOligos)):
		outFile.write(">oligo_" + str(i) + "\n")
		outFile.write(designedOligos[i] + "\n")
	outFile.close()
	return designedOligoFrequencies



if __name__ == "__main__":
	targetDirectoryPath = sys.argv[1] # path to a directory containing fasta files giving the sequences we want the oligos to hybridize to
	backgroundDirectoryPath = sys.argv[2] # path to a directory containing fasta files containing sequences we don't want the oligos to hybridize to
	targetOligoLength = int(sys.argv[3]) # desired length of oligos
	maxGapLength = int(sys.argv[4]) # desired length of oligos
	outputPath = sys.argv[5] # path to write output files
	
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

	if not os.path.exists(outputPath+"/backgroundFrequencyDict"):
		# compute background frequency dict only if it's required, otherwise it's a pretty expensive operation
		print("reading background files")
		for backgroundFile in os.listdir(backgroundDirectoryPath):
			print(backgroundFile)
			backgroundSequences, backgroundHeaders = readFastaFile(backgroundDirectoryPath + "/" + backgroundFile)
			allBackgroundSequences += backgroundSequences
			allBackgroundHeaders += backgroundHeaders
		# compute frequency of background oligos	
		backgroundFrequencyDict = countOligos(allBackgroundSequences, targetOligoLength)

		backgroundFrequencyFile = open(outputPath+"/backgroundFrequencyDict","wb")
		pickle.dump(backgroundFrequencyDict, backgroundFrequencyFile)
		backgroundFrequencyFile.close()

	else:
		backgroundFrequencyDict = pickle.load(open(outputPath + "/backgroundFrequencyDict", "rb"))
	
	designedOligoFrequencies = designOligos(allTargetSequences, backgroundFrequencyDict, targetLength, maxGapLength, outputPath)
		
	# plot background frequencies
	plt.hist(list(backgroundFrequencyDict.values()))
	plt.xlabel("Oligo Frequency")
	plt.ylabel("Frequency")
	plt.savefig(outputPath + "/backgroundFrequency.png")
	plt.close()

	# plot oligo frequency distribution vs background
	plt.hist(designedOligoFrequencies)
	plt.xlabel("Oligo Frequency")
	plt.ylabel("Frequency")
	plt.savefig(outputPath + "/designedOligoFrequency.png")
	plt.close()
	


