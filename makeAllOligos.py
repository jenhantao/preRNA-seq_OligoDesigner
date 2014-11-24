'''
Designs oligos for a pre RNA-seq selection method
'''

### imports ###
import sys
import os
import numpy as np


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

def makeOligos(targetSequences, targetLength, outputPath):
	'''
	Gives all non-unique k-mers of target length that appear in target sequences
	inputs: a list of sequences, length of k-mers, path to write output files
	outputs: writes the designed oligos to a Fasta file

	'''
	seenOligos = set()
	for i in range(len(targetSequences)):
		currentSeq = targetSequences[i]
		for j in range(len(targetSequences[i]) - targetLength):
			oligo = currentSeq[ j : j + targetLength ]
			seenOligos.add(oligo)

	# write fasta files
	oligos = list(seenOligos)
	for i in range(len(oligos)):
		outFile = open(outputPath + "/" + oligos[i] + ".fa", "w")
		for j in range(1):
			outFile.write(">" + str(j) + "\n")
			outFile.write(oligos[i] + "\n")
		outFile.close()
		
if __name__ == "__main__":
	targetDirectoryPath = sys.argv[1] # path to a directory containing fasta files giving the sequences we want the oligos to hybridize to
	targetLength = int(sys.argv[2]) # desired length of oligos
	outputPath = sys.argv[3] # path to write output files
	
	# intialize lists
	allTargetSequences = []
	allTargetHeaders = []

	# read in sequences
	print("reading target files")
	for targetFile in os.listdir(targetDirectoryPath):
		print(targetFile)	
		targetSequences, targetHeaders = readFastaFile(targetDirectoryPath + "/" + targetFile)
		allTargetSequences += targetSequences
		allTargetHeaders += targetHeaders

	print("writing oligo fasta files")
	makeOligos(targetSequences, targetLength, outputPath)
