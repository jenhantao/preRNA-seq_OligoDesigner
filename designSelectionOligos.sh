# read in inputs
targetDirectoryPath=$1 # path to a directory containing fasta files giving the sequences we want the oligos to hybridize to
targetOligoLength=$2 # desired length of oligos
maxGapLength=$3 # desired length of oligos
threshold=$4
outputPath=$5 # path to write output files


if [ ! -d $outputPath ]
then
	mkdir $outputPath
fi

# design all oligos
if [ ! -d $outputPath/oligoFastas ] # make all oligos only if the directory doesn't exist; this is an expensive step
then
	mkdir $outputPath/oligoFastas
	echo "enumerating all oligos"
	command="python makeAllOligos.py $targetDirectoryPath $targetOligoLength $outputPath/oligoFastas"
	echo $command
	$command
fi

# count background frequency using bowtie
if [ ! -d $outputPath/oligoFasta ] # create fasta files only if the directory doesn't exist; this is an expensive step
then
	export BOWTIE_INDEXES="/data/home/jenhan/analysis/oligo_design/mm10/" # TODO fix bowtie installation
	echo "performing alignments"
	for sequenceFasta in $outputPath/oligoFastas/*.fa
	do
		alignPath=${sequenceFasta%.fa}
		alignLogPath="${alignPath}_log"
		alignPath+="_align"
		command="bowtie genome -a -v 0 -f $sequenceFasta "
		echo $command
		$command > $alignPath 2> $alignLogPath
		
	done
fi


# create background frequency dict
echo "counting mapped reads"
backgroundFrequencyDictPath="$outputPath/backgroundFrequencyDict" # path to a directory containing fasta files containing sequences we don't want the oligos to hybridize to
command="python mappedReadsCounter.py $outputPath/oligoFastas $backgroundFrequencyDictPath"
echo $command
$command
exit 0



# filter oligos according to background count
echo "filtering oligos"
command="python oligoDesigner.py $targetDirectoryPath $backgroundFrequencyDictPath $mapGapLength $threshold $outputPath"
echo $command
$command

