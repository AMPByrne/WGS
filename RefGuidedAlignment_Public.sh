#!/bin/bash
set -e

#BWA Mapping of Oxford Nanopore and Illumina Reads.
#version 1.0 - changed the name from BWAMapping.sh to RefGuidedAlignment.sh after version 2.0. The usage references and script name have been updated throughout. The comments and commands have also been cleaned up as well to remove specific locations.
#version 2.0 - replaced the use of vcf2consensus.pl with genconsensus.py which uses the program iVar to generate the final consensus. Dependencies updated accordingly. By default, genconsensus.py will call a degenerate base at a MinFreq of 0.8, which is equivalent to the bases that make up 80% of the depth at a position. E.g. if you have 10 reads at a position, consisting of 7 As and 3 Gs, genconsensus.py will call an R as 80% of the reads consist of As AND Gs. See iVar manual for more information: https://andersen-lab.github.io/ivar/html/manualpage.html. Also by default, any positions without any read coverage, or read coverage not of a suitable quality, will be replaced by dashes (-) in the final consensus.
#version 2.1 - this script has been altered to be used with the long reads that are produced by Oxford Nanopore reads. 
#version 3.0 - integrated the ability to process Oxford Nanopore and Illumina reads into one script through options.Changed the way the script is used, by integrating options. No longer have to provide a sample name, this is taken from the readfile name. The output location has been changed so that unless stated otherwise, it will just output the data to where the read files are stored. Just generally rearranged the script as well to make more sense. Added the ability for users to alter the variables for genconsensus.py. Note - in version 2.0 of genconsensus.py, the default minimum frequency threshold has been changed from 0.8 to 0, to ensure the base that makes up the majority is used to call the consensus.
#
#
#Required programs/tools:
#BWA
#Samtools
#bcftools
#iVar
#genconsensus.py

VERSION=3
#time the process
Start=$(date +%s)
RunDate=$(date +%d-%m-%y)

#Default options
Cores=$(grep -c ^processor /proc/cpuinfo)
#iVar options
MinFreq=0
MinDepth=1
NoCoverage=-

#Define a function for usage
usage()
{ echo "
RefGuidedAlignment.sh Version $VERSION
Usage: bash $0 [options] <read(s).fastq>
	-r path to reference genome file [Required]
	-s sequencer used, Illumina or ONT [Required]
	-o path to output location [optional] default = location of read file(s)
	-c cores [optional] default = all
	-t the minimum frequency at which a nucleotide will be used to call the consensus. INT 0 - 1.[optional] deault = $MinFreq
	-m minimum read depth at which iVar call a consensus [optional] default = $MinDepth
	-n (N or -) character to print if less than minimum coverage (-m) [optional] default = $NoCoverage

"
}
if [ "$1" = help ] || [ -z "$1" ]; then
usage
	exit
fi
#Parse the options
while getopts 'r:s:o:c:t:m:n' opt ; do
	case $opt in
		r) Reference=$OPTARG ;;
		s) Sequencer=$OPTARG ;;
		o) OutputLocation=$OPTARG ;;
		c) Cores=$OPTARG ;;
		t) MinFreq=$OPTARG ;;
		m) MinDepth=$OPTARG ;;
		n) NoCoverage=$OPTARG ;;

	esac
done
#skip over the processed options
shift $((OPTIND-1))

#check for mandatory positional parameters
if [[ "$Reference" == "" ]]; then
usage && echo "*****Error - there needs to be a reference sequence*****

"
	exit
fi
if [[ "$Sequencer" == "" ]]; then
usage && echo "*****Error - you need to state whether you used an Illumina or ONT sequencer*****

"
	exit
fi

#Check for read files and how many have been provided and if this matches the Sequencer provided
if [ ! -s "$1" ] ; then
usage && echo "*****Error - you need to provide at least one read file*****

"
	exit
fi

if [ "$1" ] && [ "$2" ] ; then
ReadType=Illumina
R1reads="$1"
R2reads="$2"
fi

if [ "$1" ] && [ ! -s "$2" ] ; then
ReadType=ONT
R1reads="$1"
fi

if [ "$ReadType" == "Illumina" ] && [ "$Sequencer" == "ONT" ] ; then
usage && echo "*****Error - you have selected the Sequencer ONT but provided 2 read files - please check Sequencer*****

"
	exit
fi

if [ "$ReadType" == "ONT" ] && [ "$Sequencer" == "Illumina" ] ; then
usage && echo "*****Error - you have selected the Sequencer Illumina but provided 1 read file - please check Sequencer*****

"
	exit
fi

echo "Pre-processing checks OK, let's begin processing"

#Get directory/sample/reference names from file names
Filename=$(basename "$1")
Samplename=${Filename%%_*}

Ref=$(basename "$Reference")

if [ "$OutputLocation" == "" ] ; then
	DIR=$(dirname "$1")
else
	DIR="$OutputLocation"
fi

OutputDir="$DIR"/"$Samplename"_RefGuidedAlignment_"$Ref"_"$RunDate"

#See if OutputDir already exists
test -d $OutputDir && usage && echo "*****Error - Output directory already exists, you need to fix this*****" && exit || echo Output Directory being created
mkdir "$OutputDir" && cd "$OutputDir"
cp "$Reference" "$OutputDir"
Reference="$OutputDir"/"$Ref"
#
#Map reads to reference and output a sam file
echo Aligning reads to reference sequence
if [ "$Sequencer" == "Illumina" ] ; then 
echo "Mapping reads to $Ref using bwa mem"
bwa index "$Reference"
bwa mem -t"$Cores" "$Reference" "$R1reads" "$R2reads" | samtools fixmate -O bam - "$Samplename".bam
fi

if [ "$Sequencer" == "ONT" ] ; then
echo "Mapping reads to $Ref using minimap2"
minimap2 -ax map-ont -t "$Cores" "$Reference" "$R1reads" | samtools fixmate -O bam - "$Samplename".bam
fi

#sort the bam file
echo Sorting bam file
samtools sort "$Samplename".bam -o "$Samplename".sorted.bam

#generate mapping stats
echo Generating mapping stats
samtools flagstat "$Samplename".sorted.bam > "$Samplename".sorted.bam.STATS.pre-rmdup.txt

#remove duplicate reads
echo Removing duplicate reads
samtools rmdup -sS "$Samplename".sorted.bam "$Samplename".sorted.rmdup.bam

#generate mapping stats post duplicate removal
echo Generating mapping stats post duplicate removal
samtools flagstat "$Samplename".sorted.rmdup.bam > "$Samplename".sorted.bam.STATS.rmdup.txt

#index the bam file
echo Indexing bam file
samtools index -b "$Samplename".sorted.rmdup.bam
#
#create a consensus sequence from the rmdup.sorted.bam file using genconsensus.py
echo Generating consesus sequence
mkdir genconsensus_Results
cd genconsensus_Results
#change the below to provide the location of genconsensus.py
python /path/to/genconsensus.py -t "$MinFreq" -m "$MinDepth" -n "$NoCoverage" -r "$Reference" -b ../"$Samplename".sorted.rmdup.bam
mv final_consensus.fasta ../"$Samplename".fasta
cd ..		
rm -r genconsensus_Results
#
#just clearing up now and move the final files to the output location
echo Time to clean up
rm "$Samplename".bam
rm "$Samplename".sorted.bam
rm "$Reference".*

End=$(date +%s)
TimeTaken=$((End-Start))
echo "All finished,you're results are here:' $OutputDir"
echo | awk -v D=$TimeTaken '{printf "RefGuidedAlignment.sh took %02d'h':%02d'm':%02d's'\n",D/(60*60),D%(60*60)/60,D%60}'
