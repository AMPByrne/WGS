#!/bin/bash
set -e

#version 0 - The initial version of this script was based on FluSeqID (https://github.com/ellisrichardj/FluSeqID) but was re-written by James Seekings to run on CentOS rather than Biolinux.
#verion 1.0 - Velvet was replaced with SPAdes due to ongoing issues. This script was then renamed to denovoAssembly.sh as it had changed significantly form the original FluSeqID script.
#version 2.0 - Removed the use vcf2consensus.pl for generating the final consensus. Instead this has now been replaced with genconsensus.py, which uses the program iVar to generate the final consensus. Dependencies updated accordingly. Fixed the inclusion of the data in the results directory.
#version 3.0 - Updated to work with genconsensus.py version 2.0, which allows the user to define variables for iVar. The usage function has been updated accordingly. The flag for threads has been changed from -t to -T to avoid conflict with iVar flags. Note - in version 2.0 of genconsensus.py, the default minimum frequency threshold has been changed from 0.8 to 0, to ensure the base that is observed in the majority of reads (accounting for read quality) is used to call the consensus. The default minimum read depth is still 1, and the default character used for nucleotides without coverage is still -.
#Version 4.0 - changed syntax for subsampling from [$a -gt $b] to (($a > $$b)) this should fix subsampling issues (TO BE REVIEWED)


#This script does the following:
#Remove host genome from raw illumina reads
#Subsample the bam if necesary
#Output non host paired reads from the bam file
#Denovo assemble with SPAdes
#Blast the contigs against a local database
#Use the top hit(s) from blast to produse a reference genome
#check that there were any blast hits
#Map the raw illumina reads agianst the reference genoome
#Produce a new consensus and use this as a reference for subsequent mapping

#Dependencies:
#BWA
#samtools
#bamtools
#SPAdes
#blast+
#picard
#GATK
#iVar

VERSION=3.0
#time the process
Start=$(date +%s)
begin=$(date '+%x %R')
RunDate=$(date +%d-%m-%y)

#Defaults for options
PathToHOST=/path/to/GCF_000003025.6_Sscrofa11.1_genomic.fna.gz
PathToSearchData=/path/to/Viral_Reference_Databases
Blast_e_value=0.0001
threads=$(grep -c ^processor /proc/cpuinfo)
MinFreq=0
MinDepth=1
NoCoverage=-


#define a function for usage
usage()
{ echo "
denovoAssembly.sh Version $VERSION
Usage: bash $0 [options] -G <HostGenome> <Read1.fastq.gz> <Read2.fastq.gz>
	help Displays this usage
	-G path to HOST genome File [Default] 
	-r path to Reference Database Directory [optional] default = $PathToSearchData
	-e Blast e value [optional] default = $Blast_e_value
	-T number of threads/cores to use [optional] default = all
	-t the minimum frequency at which a nucleotide will be used to call the consensus. INT 0 - 1.[optional] deault = $MinFreq
	-m minimum read depth at which iVar call a consensus [optional] default = $MinDepth
	-n (N or -) character to print if less than minimum coverage (-m) [optional] default = $NoCoverage
";
}

# parse the options
while getopts 'G:e:r:t:' opt ; do
  case $opt in
    G) PathToHOST=$OPTARG ;;
    e) Blast_e_value=$OPTARG ;;
    r) PathToSearchData=$OPTARG ;;
    T) threads=$OPTARG ;;
    t) MinFreq=$OPTARG ;;
    m) MinDepth=$OPTARG ;;
    n) NoCoverage=$OPTARG ;;
    
  esac
done
# skip over the processed options
shift $((OPTIND-1)) 


# check for mandatory positional parameters - Read Files
if [[ $1 = help ]]; then
usage
  exit
fi
if [ $# -ne 2 ]; then
usage && echo "*****Error - There needs to be 2 read files*****"
  exit
fi

LEFT="$(readlink -f "$1")"
RIGHT="$(readlink -f "$2")"

#Check read files exist
if [ ! -s "$LEFT" ] ; then
usage && echo "*****Error - reads1 has a problem, you need to investigate*****" && exit
fi

if [ ! -s "$RIGHT" ] ; then
usage && echo "*****Error - reads2 has a problem, you need to investigate*****" && exit
fi

# Get directory/sample/reference/host names from file names
	sfile1=$(basename "$LEFT")
	sfile2=$(basename "$RIGHT")
	samplename=${sfile1%%_*}

	ref=$(basename "$PathToSearchData")

	host=$(basename "$PathToHOST")
	hostname=${host%%.*}
#create the output folder in same location as data
	DIR=$(dirname "$LEFT")


# Create output OutputDirectory
OutputDir="$DIR"/"$samplename"_SPAdes_"$ref"_"$RunDate"

#See if OutputDir already exists
test -d $OutputDir && usage && echo *****Error - OutputDirectory already exists, you need to fix this***** && exit || echo OutputDirectory being created
mkdir "$OutputDir"

#can add version info at this point. BM removed as not needed at the moment

echo "Getting host genome"

#Generate bwa index of host genome if it doesn't already exist
if [ -s "$PathToHOST".sa ] && [ -s "$PathToHOST".amb ] && [ -s "$PathToHOST".ann ] && [ -s "$PathToHOST".bwt ] && [ -s "$PathToHOST".pac ]
	then
		echo "Using existing host index"
	else
		echo "Generating local index of host genome"
	ln -s "$(readlink -f "$PathToHOST")" "$OutputDir"/"$hostname".fa
	bwa index "$OutputDir"/"$hostname".fa
	PathToHOST="$OutputDir"/"$hostname".fa
fi

echo "Mapping raw reads to Host"

#Map reads1 and reads2 to the host genome. Will output a bam file
bwa mem -t "$threads" "$PathToHOST" "$LEFT" "$RIGHT" | samtools view -@ "$threads" -b -f 4 -o "$OutputDir"/"$samplename"_nonHost.bam -;

#Subsample bam if necessary
NonHostReads=$(samtools view -c "$OutputDir"/"$samplename"_nonHost.bam)
MaxReads=2000000
if (( $NonHostReads > $Maxreads ))
	then	Sub=$(echo "scale=2; 321+$MaxReads/$NonHostReads" | bc)
			samtools view -@ "$threads" -bs "$Sub" "$OutputDir"/"$samplename"_nonHost.bam -o "$OutputDir"/"$samplename"_nonHost_subsample.bam
		SPAdesInBam="$OutputDir"/"$samplename"_nonHost_subsample.bam
	else	SPAdesInBam="$OutputDir"/"$samplename"_nonHost.bam
fi

#Extract all paired reads that are NOT mapped to the host genome
samtools sort -@ "$threads" -n -O BAM "$SPAdesInBam" | samtools fastq -@ "$threads" -1 "$OutputDir"/"$samplename"_NotHostReads1.fastq -2 "$OutputDir"/"$samplename"_NotHostReads2.fastq -s "$OutputDir"/"$samplename"_NotHostReadsSingletons.fastq -;

echo "Host Removed"

echo "NonHostReads="$NonHostReads""

echo "Starting SPAdes Assembler"

#Start SPAdes denovo assembly of non host readsmake
#SPAdes.py can be changed to SPAdes if you create symbolic link
spades.py --threads "$threads" --only-assembler --pe1-1 "$OutputDir"/"$samplename"_NotHostReads1.fastq --pe1-2 "$OutputDir"/"$samplename"_NotHostReads2.fastq -o "$OutputDir"/SPAdes/

#Copy the contig.fasta
cp "$OutputDir"/SPAdes/contigs.fasta "$OutputDir"

#Run Blast of output contigs against database
cd "$OutputDir"
echo "BLASTing contigs against selected database(s)"
for file in "$PathToSearchData"/*.fasta
do
	seg=$(basename "$file")
	segname=${seg%%.*}

	blastn -db "$file" -query contigs.fasta -out "$samplename"_"$segname"_crunch.txt -evalue "$Blast_e_value" -max_target_seqs 5 -outfmt 6 -num_threads "$threads" &

done
wait

#Extract the single sequence  from blast output corresponding to the highest scoring match for each genome segment

for crunch in *crunch.txt
do
	topseg=${crunch%%.*}
	part=${topseg#*_}
	segname=${part%_*}

	sort -k12,12 -rn "$crunch" | head -1 - | awk '{print $2}' - > "$topseg"_match.txt

		blastdbcmd -db "$PathToSearchData"/"$segname".fasta -entry_batch "$topseg"_match.txt > "$segname"_match.fas &
done
wait
cat *_match.fas > top_matches.fa

#Check file size of top_matches.fa to make sure it has data
if [ ! -s top_matches.fa ] ; then
echo "*****Error - top_matches.fa is empty. No contigs had blast matches*****" && exit
fi

#Copy all the reference names to add to file  later
REFS=$( grep '>' top_matches.fa )

	#Map raw data to selected reference sequences and generate new consensus (iteratively)
	rfile="$(readlink -f top_matches.fa)"

		ref=$(basename "$rfile")
		refname=${ref%%_*}
		reffile=${ref%%.*}

	iter=4
	count=1

	mkdir "$samplename"_IterMap"$iter"
	cd "$samplename"_IterMap"$iter"

#Log commands that are run
echo -e "$now\n\denovovassembly.sh running with $threads cores\n\The following commands were run: \n" > "$samplename"_IterMap"$iter".log

# Define function to log commands and then run them
	LogRun(){
	echo -e "\n$@" >> "$samplename"_IterMap"$iter".log
	eval "$@"
	}

		while (($count <= $iter))
		do
		# Set mapping stringency for first iteration - could reduce it
		if [ $count == 1 ] && [ $iter != 1 ]; then
			mem=19
			mmpen=4
			gappen=6
		else # bwa mem defaults
			mem=19
			mmpen=4
			gappen=6
		fi

		# mapping to original reference or most recently generated consensus
		LogRun bwa index "$rfile"
		LogRun samtools faidx "$rfile"
		#Trying to multithread java 'java -XX:ParallelGCThreads=<num of  threads> -jar'
		#BM changed so picard works without sending for jar file
		LogRun picard -XX:ParallelGCThreads="$threads" CreateSequenceDictionary R="$rfile" O=${rfile%%.*}.dict
#09 Oct 2018 added F4 to try and make samtools view output only mapped reads in each iteration - took out F4 again as stats no longer said what % of reads were virus
		LogRun bwa mem -T10 -t "$threads" -k "$mem" -B "$mmpen" -O "$gappen" -R '"@RG\tID:"$samplename"\tSM:"$samplename"\tLB:"$samplename""' "$rfile" "$LEFT" "$RIGHT" |
		samtools view -@ "$threads" -u - | samtools sort -@ "$threads" -o "$samplename"_"$reffile"-iter"$count"_map_sorted.bam
		samtools index -@ "$threads" "$samplename"_"$reffile"-iter"$count"_map_sorted.bam
		#LogRun java -XX:ParallelGCThreads="$threads" -jar /home/$USER/mnt/VI6Bioinformatics/Central_Pipelines/denovoAssembly/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt "$threads" -R "$rfile" -I "$samplename"_"$reffile"-iter"$count"_map_sorted.bam -o indel"$count".list
		#Trying to use the gatk.py script rather than the .jar file to see if that solves the issues with accessing the .jar file
		LogRun python '/path/to/gatk.py' -T RealignerTargetCreator -nt "$threads" -R "$rfile" -I "$samplename"_"$reffile"-iter"$count"_map_sorted.bam -o indel"$count".list
		#LogRun java -XX:ParallelGCThreads="$threads" -jar /home/$USER/mnt/VI6Bioinformatics/Central_Pipelines/denovoAssembly/GenomeAnalysisTK.jar -T IndelRealigner -R "$rfile" -I "$samplename"_"$reffile"-iter"$count"_map_sorted.bam -targetIntervals indel"$count".list -maxReads 50000 -o "$samplename"_"$reffile"-iter"$count"_realign.bam
		LogRun python '/path/to/gatk.py' -T IndelRealigner -R "$rfile" -I "$samplename"_"$reffile"-iter"$count"_map_sorted.bam -targetIntervals indel"$count".list -maxReads 50000 -o "$samplename"_"$reffile"-iter"$count"_realign.bam
		rm "$samplename"_"$reffile"-iter"$count"_map_sorted.bam
		rm "$samplename"_"$reffile"-iter"$count"_map_sorted.bam.bai

	if [ $count == $iter ]; then
		#generate and correctly label consensus using cleaned bam on final iteration, removing duplicates
		#using samtools markdup to remove duplicates
		LogRun samtools sort -n -@ "$threads" "$samplename"_"$reffile"-iter"$count"_realign.bam -o - | \
		samtools fixmate -r -m 	-@ "$threads" - - | \
		samtools sort -@ "$threads" -o - | \
		samtools markdup -r -s -O BAM -@ "$threads" - "$samplename"_"$reffile"-iter"$count"_clean_mapOnly.bam ;
		
		samtools index -@ "$threads" "$samplename"_"$reffile"-iter"$count"_clean_mapOnly.bam

        #replaced vcf2consensus.pl with genconsensus.py here
        mkdir genconsensus_Results
		cd genconsensus_Results
		python /path/to/genconsensus.py  -t "$MinFreq" -m "$MinDepth" -n "$NoCoverage" -r ../"$rfile" -b ../"$samplename"_"$reffile"-iter"$count"_clean_mapOnly.bam
        	LogRun mv final_consensus.fasta ../"$samplename"_"$reffile"-iter"$count"_consensus.fasta
		cd ..		
		rm -r genconsensus_Results
              
		# generate mapping statistics
		LogRun samtools flagstat "$samplename"_"$reffile"-iter"$count"_realign.bam > "$samplename"_"$reffile"-iter"$count"_MappingStats.txt
		

	else
        #kept vcf2consensus.pl here. This is just generating a reference to be used in the mapping for the next iteration. genconsensus.py will remove any sequence without reads, which could negatively impact the next mapping iteration.
        	LogRun bcftools mpileup -L 10000 -Q1 -AEpf "$rfile" "$samplename"_"$reffile"-iter"$count"_realign.bam | bcftools call --threads "$threads" -c - > "$samplename"_"$reffile"-iter"$count".vcf
		LogRun perl /path/to/vcf2consensus.pl  consensus -f "$rfile" "$samplename"_"$reffile"-iter"$count".vcf |
		sed '/^>/ s/-iter[0-9]//;/^>/ s/$/'-iter"$count"'/' - > "$samplename"_"$reffile"-iter"$count"_consensus.fasta

		# mapping statistics
		LogRun samtools flagstat "$samplename"_"$reffile"-iter"$count"_realign.bam > "$samplename"_"$reffile"-iter"$count"_MappingStats.txt
		rfile="$samplename"_"$reffile"-iter"$count"_consensus.fasta

	fi

		((count=count+1))
		echo "New Consensus: "$rfile""
	done


	complete=$(date '+%x %R')
	LogRun echo "Started at "$begin""
	LogRun echo "Completed at "$complete""
LogRun echo "Results are in: "$OutputDir"/"$samplename"_FinalResults"



	echo -e "Completed processing "$complete""  >> "$samplename"_IterMap"$iter".log

	# Clean up intermediate files
	cd ..
	#Define Final Results folder with a date
	FinalResults="$samplename"_FinalResults_"$RunDate"
	mkdir "$FinalResults"
	mv "$samplename"_IterMap"$iter"/"$samplename"_"$reffile"-iter4_MappingStats.txt  "$FinalResults"/"$samplename"_"$reffile"-iter4_MappingStats.txt
	mv "$samplename"_IterMap"$iter"/"$samplename"_"$reffile"-iter4_consensus.fasta "$FinalResults"/"$samplename"_"$reffile"-iter4_consensus.fasta
	mv "$samplename"_IterMap"$iter"/"$samplename"_"$reffile"-iter3_consensus.fasta "$FinalResults"/"$samplename"_"$reffile"-iter3_consensus.fasta
	mv "$samplename"_IterMap"$iter"/"$samplename"_"$reffile"-iter4_clean_mapOnly.bam "$FinalResults"/"$samplename"_"$reffile"-iter4_clean_mapOnly.bam
	mv "$samplename"_IterMap"$iter"/"$samplename"_"$reffile"-iter4_clean_mapOnly.bam.bai "$FinalResults"/"$samplename"_"$reffile"-iter4_clean_mapOnly.bam.bai
	mv "$samplename"_IterMap"$iter"/"$samplename"_IterMap"$iter".log "$FinalResults"/"$samplename"_IterMap"$iter".log

End=$(date +%s)
TimeTaken=$((End-Start))
echo  | awk -v D=$TimeTaken '{printf "Performed denovoassembly Analysis in: %02d'h':%02d'm':%02d's'\n",D/(60*60),D%(60*60)/60,D%60}'
