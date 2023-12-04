#!/bin/bash
set -e
#Alex Byrne
#version 1.0	30Jun2023	The initial version of this script to perform hybrid denovo/reference-guided assembly on nanopore sequencing reads. This pipeline used denovoAssembly.sh version 5.0 as a basis and was then adapted to be compatible with long-read nanopore data - as such certain tools/steps were changed or omitted:
	# -BWA-MEM was swapped for Minimap2 for read alignment
	# -The SPAdes step was removed as nanopore reads are of sufficient length to be 	appropriate for querying BLAST databases. However, this does create large BLAST hit files 	(MBs versus kbs that were obtained from short-read pipeline. This may be optimised in the	future through sub-sampling of the non-host reads but at the moment this will work.	
	# -vcf2consensus.pl and iVar were swapped entirely for Medaka for generating the intermediate and final consensus sequences. Medaka is designed specifically for calling consensus sequences from nanopore data and is therefore most appropriate for this purpose. 

#This script does the following:
#Aligns reads to the ONT DNA control (lambda phage DNA) sequence used as an internal control in V14 kits and determines the read depth and coverage - this is useful to troubleshoot failed runs
#Remove host genome from raw nanopore reads using Minimap2
#Output non-host reads from the bam file using Samtools
#BLAST the non-host reads against a local database using blastn
#Use the top hit(s) from blast to produce to select a reference genome
#Map the raw nanopore reads agianst the reference genoome using Minimap2 and Samtools
#Produce a consensus from this alignment using Medaka and BCFTools which is then used for iterative mapping of the raw reads.
#For the final consensus sequence generation, regions with read depths less than the defined minimum read coverage (currently 10 reads) are masked.

#It's worth noting that the removal of the SPAdes assembly step makes the script much more computationally efficient meaning it can be comfortably run on a less powerful machine - we've had success running this pipeline with only 4 CPU, 16 GB RAM, whereas the denovoAssembly.sh pipline requires 16 CPU, 64 GB RAM.

#Dependencies:
#Minimap2
#samtools
#bcftools
#blast+
#medaka

VERSION=1.0
#time the process
Start=$(date +%s)
begin=$(date '+%x %R')
RunDate=$(date +%d-%m-%y)

#Defaults for options
PathToHOST=/path/to/GCF_000003025.6_Sscrofa11.1_genomic.fna.gz
PathToSearchData=/path/to/Viral_Reference_Databases/
MedakaModel=r1041_e82_400bps_hac_variant_v4.2.0 
Blast_e_value=0.0001
threads=$(grep -c ^processor /proc/cpuinfo)
MinDepth=5
PathToDCS=/path/to/ONT_DCS_Sequence.fasta


#define a function for usage
usage()
{ echo "
denovoAssembly_ONT.sh Version $VERSION
Usage: bash $0 [options] -G <HostGenome> -r <ReferenceDatabase> <Reads.fastq.gz> 
	help Displays this usage
	-G path to HOST genome File [Default] 
	-r path to Reference Database Directory [optional] default = $PathToSearchData
	-e Blast e value [optional] default = $Blast_e_value
	-T number of threads/cores to use [optional] default = all
	-m minimum read depth at which a consensus sequence is called. Region with read coverage below this will be masked [optional] default = $MinDepth
	"
}

# parse the options
while getopts 'G:e:r:t:' opt ; do
  case $opt in
    G) PathToHOST=$OPTARG ;;
    e) Blast_e_value=$OPTARG ;;
    r) PathToSearchData=$OPTARG ;;
    T) threads=$OPTARG ;;
    m) MinDepth=$OPTARG ;;
    
  esac
done
# skip over the processed options
shift $((OPTIND-1)) 


# check for mandatory positional parameters - Read Files
if [[ $1 = help ]]; then
usage
  exit
fi
if [ $# -ne 1 ]; then
usage && echo "*****Error - There needs to be 1 read files*****"
  exit
fi

LEFT="$(readlink -f "$1")"

#Check read files exist
if [ ! -s "$LEFT" ] ; then
usage && echo "*****Error - reads have a problem, you need to investigate*****" && exit
fi

# Get directory/sample/reference/host names from file names
	sfile1=$(basename "$LEFT")
	samplename=${sfile1%%_*}

	ref=$(basename "$PathToSearchData")

	host=$(basename "$PathToHOST")
	hostname=${host%%.*}
#create the output folder in same location as data
	DIR=$(dirname "$LEFT")


# Create output OutputDirectory
OutputDir="$DIR"/"$samplename"_ONT_"$ref"_"$RunDate"

#See if OutputDir already exists
test -d $OutputDir && usage && echo *****Error - OutputDirectory already exists, you need to fix this***** && exit || echo OutputDirectory being created
mkdir "$OutputDir"

echo "Mapping raw reads to DNA Control Sample Sequence"
cp "$PathToDCS" "$OutputDir"/

minimap2 -ax map-ont -s 10 -t "$threads" "$PathToDCS" "$LEFT" | samtools view -@ "$threads" -u - | samtools sort -@ "$threads" -o "$OutputDir"/"$samplename"_DCS_map_sorted.bam

samtools coverage "$OutputDir"/"$samplename"_DCS_map_sorted.bam > "$OutputDir"/"$samplename"_DCS_map_sorted.coverage

samtools depth "$OutputDir"/"$samplename"_DCS_map_sorted.bam > "$OutputDir"/"$samplename"_DCS_map_sorted.depth

mkdir "$OutputDir"/DCS_Outputs

mv "$OutputDir"/"$samplename"_DCS_map_sorted.bam "$OutputDir"/DCS_Outputs/"$samplename"_DCS_map_sorted.bam
mv "$OutputDir"/"$samplename"_DCS_map_sorted.coverage "$OutputDir"/DCS_Outputs/"$samplename"_DCS_map_sorted.coverage
mv "$OutputDir"/"$samplename"_DCS_map_sorted.depth "$OutputDir"/DCS_Outputs/"$samplename"_DCS_map_sorted.depth

echo "Mapping raw reads to Host"

#Map reads to the host genome. Will output a bam file
minimap2 -ax map-ont -t "$threads" "$PathToHOST" "$LEFT" | samtools view -@ "$threads" -b -f 4 -o "$OutputDir"/"$samplename"_nonHost.bam -;

#Extract all reads that are NOT mapped to the host genome
samtools sort -@ "$threads" -n -O BAM "$OutputDir"/"$samplename"_nonHost.bam | samtools fasta -@ "$threads" -0 "$OutputDir"/"$samplename"_NotHostReads.fasta -;


echo "Host Removed"

#Run Blast of output contigs against database
echo "BLASTing contigs against selected database(s)"
for file in "$PathToSearchData"/*.fasta
do
	seg=$(basename "$file")
	segname=${seg%%.*}

	blastn -db "$file" -query "$OutputDir"/"$samplename"_NotHostReads.fasta -out "$OutputDir"/"$samplename"_"$segname"_crunch.txt -evalue "$Blast_e_value" -max_target_seqs 5 -outfmt 6 -num_threads "$threads" &

done
wait

#Extract the single sequence from blast output corresponding to the highest scoring match for each genome segment
cd "$OutputDir"
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
cd

#Check file size of top_matches.fa to make sure it has data
if [ ! -s "$OutputDir"/top_matches.fa ] ; then
echo "*****Error - top_matches.fa is empty. No contigs had blast matches*****" && exit
fi

#Removing / from headers as this breaks medaka
sed -i 's+/+_+g' "$OutputDir"/top_matches.fa

#Remove any degenerate bases from the reference sequences as this also breaks medaka - replace all with N as medaka can de-ambiguate this but not other degenerate bases
sed -i '/^[^>]/s/[RYSWKMBDHV]/N/g' "$OutputDir"/top_matches.fa

#Copy all the reference names to add to file later
REFS=$( grep '>' "$OutputDir"/top_matches.fa )

	#Map raw data to selected reference sequences and generate new consensus (iteratively)
	rfile="$(readlink -f "$OutputDir"/top_matches.fa)"

		ref=$(basename "$rfile")
		refname=${ref%%_*}
		reffile=${ref%%.*}

	iter=4
	count=1

	mkdir "$OutputDir"/"$samplename"_IterMap"$iter"
	cd "$OutputDir"/"$samplename"_IterMap"$iter"
	qt="'"
	string='ALT="."'
	

# Define function to log commands and then run them
	LogRun(){
	echo -e "\n$@" >> "$samplename"_IterMap"$iter".log
	eval "$@"
	}

		while (($count <= $iter))
		do
		# mapping to original reference or most recently generated consensus		
		LogRun minimap2 -ax map-ont -s 10 -t "$threads" "$rfile" "$LEFT" | samtools view -@ "$threads" -u - | samtools sort -@ "$threads" -o "$samplename"_"$reffile"-iter"$count"_map_sorted.bam
		LogRun samtools index -@ "$threads" "$samplename"_"$reffile"-iter"$count"_map_sorted.bam
				

	if [ $count == $iter ]; then
		#generate and correctly label consensus using cleaned bam on final iteration. Adapted from the Epi2Me wf-flu pipeline (https://github.com/epi2me-labs/wf-flu/blob/master/main.nf). For the final consensus this will mask any positions with depth less than the MinDepth variable defined above (currently 10). For all consensus calling for ONT the recommended tool is Medaka (https://github.com/nanoporetech/medaka/tree/master).
		LogRun medaka consensus "$samplename"_"$reffile"-iter"$count"_map_sorted.bam "$samplename"_"$reffile"-iter"$count"_map_sorted.hdf --model "$MedakaModel"
		LogRun medaka variant --gvcf "$rfile" "$samplename"_"$reffile"-iter"$count"_map_sorted.hdf "$samplename"_"$reffile"-iter"$count"_map_sorted.vcf --verbose --ambig_ref
		LogRun medaka tools annotate --debug --pad 25 "$samplename"_"$reffile"-iter"$count"_map_sorted.vcf "$rfile" "$samplename"_"$reffile"-iter"$count"_map_sorted.bam "$samplename"_"$reffile"-iter"$count"_map_sorted.annotated.vcf
		LogRun bcftools filter -e "$qt""$string""$qt" "$samplename"_"$reffile"-iter"$count"_map_sorted.annotated.vcf | bcftools filter -o "$samplename"_"$reffile"-iter"$count"_map_sorted.annotated.filtered.vcf -O v -e "INFO/DP<"$MinDepth"" -
		LogRun samtools depth -aa "$samplename"_"$reffile"-iter"$count"_map_sorted.bam -Q 20 -q 20 > "$samplename"_"$reffile"-iter"$count"_map_sorted.depth
		LogRun "awk -v depth=\"$MinDepth\" '{if (\$3<depth) print \$1\"\t\"\$2}' "$samplename"_"$reffile"-iter"$count"_map_sorted.depth > "$samplename"_"$reffile"-iter"$count"_map_sorted.depth.mask"
		LogRun bgzip "$samplename"_"$reffile"-iter"$count"_map_sorted.annotated.filtered.vcf
		LogRun tabix "$samplename"_"$reffile"-iter"$count"_map_sorted.annotated.filtered.vcf.gz
		LogRun bcftools consensus --mask "$samplename"_"$reffile"-iter"$count"_map_sorted.depth.mask --mask-with 'n' --mark-del '-' --mark-ins lc --fasta-ref "$rfile" -o "$samplename"_"$reffile"-iter"$count"_consensus.fasta "$samplename"_"$reffile"-iter"$count"_map_sorted.annotated.filtered.vcf.gz
		LogRun samtools flagstat "$samplename"_"$reffile"-iter"$count"_map_sorted.bam > "$samplename"_"$reffile"-iter"$count"_MappingStats.txt
		
	else
        #generate and correctly label consensus using cleaned bam on intermediate iterations. Adapted from the Epi2Me wf-flu pipeline (https://github.com/epi2me-labs/wf-flu/blob/master/main.nf). For these interediate consensus sequences this won't mask any low coverage regions as otherwise these will be maintained in the next iteration and cause improer mapping. 
		LogRun medaka consensus "$samplename"_"$reffile"-iter"$count"_map_sorted.bam "$samplename"_"$reffile"-iter"$count"_map_sorted.hdf --model "$MedakaModel"
		LogRun medaka variant --gvcf "$rfile" "$samplename"_"$reffile"-iter"$count"_map_sorted.hdf "$samplename"_"$reffile"-iter"$count"_map_sorted.vcf --verbose --ambig_ref
		LogRun medaka tools annotate --debug --pad 25 "$samplename"_"$reffile"-iter"$count"_map_sorted.vcf "$rfile" "$samplename"_"$reffile"-iter"$count"_map_sorted.bam "$samplename"_"$reffile"-iter"$count"_map_sorted.annotated.vcf
		LogRun bcftools filter -e "$qt""$string""$qt" "$samplename"_"$reffile"-iter"$count"_map_sorted.annotated.vcf | bcftools filter -o "$samplename"_"$reffile"-iter"$count"_map_sorted.annotated.filtered.vcf -O v -e "INFO/DP<"$MinDepth"" -
		LogRun bgzip "$samplename"_"$reffile"-iter"$count"_map_sorted.annotated.filtered.vcf
		LogRun tabix "$samplename"_"$reffile"-iter"$count"_map_sorted.annotated.filtered.vcf.gz
		LogRun bcftools consensus --fasta-ref "$rfile" -o "$samplename"_"$reffile"-iter"$count"_consensus.fasta "$samplename"_"$reffile"-iter"$count"_map_sorted.annotated.filtered.vcf.gz
		LogRun samtools flagstat "$samplename"_"$reffile"-iter"$count"_map_sorted.bam > "$samplename"_"$reffile"-iter"$count"_MappingStats.txt
		LogRun rfile="$samplename"_"$reffile"-iter"$count"_consensus.fasta
       
	fi

		((count=count+1))
		echo "New Consensus: "$rfile""
		done
	
#Log commands that are run
echo -e "$now\n\denovovassembly.sh running with $threads cores\n\The following commands were run: \n" > "$samplename"_IterMap"$iter".log


	complete=$(date '+%x %R')
	LogRun echo "Started at "$begin""
	LogRun echo "Completed at "$complete""
LogRun echo "Results are in: "$OutputDir"/"$samplename"_FinalResults"



	echo -e "Completed processing "$complete""  >> "$samplename"_IterMap"$iter".log

	#Define Final Results folder with a date and move the files to be kept to it.
	databaseID=$(basename "$PathToSearchData" | cut -d _ -f 1 -)
	FinalResults="$samplename"_FinalResults_"$databaseID"_Database_"$RunDate"
	mkdir "$OutputDir"/"$FinalResults"
	mv "$OutputDir"/"$samplename"_IterMap"$iter"/"$samplename"_"$reffile"-iter4_MappingStats.txt  "$OutputDir"/"$FinalResults"/"$samplename"_"$reffile"-iter4_MappingStats.txt
	mv "$OutputDir"/"$samplename"_IterMap"$iter"/"$samplename"_"$reffile"-iter4_consensus.fasta "$OutputDir"/"$FinalResults"/"$samplename"_"$reffile"-iter4_consensus.fasta
	mv "$OutputDir"/"$samplename"_IterMap"$iter"/"$samplename"_"$reffile"-iter3_consensus.fasta "$OutputDir"/"$FinalResults"/"$samplename"_"$reffile"-iter3_consensus.fasta
	mv "$OutputDir"/"$samplename"_IterMap"$iter"/"$samplename"_"$reffile"-iter4_map_sorted.bam "$OutputDir"/"$FinalResults"/"$samplename"_"$reffile"-iter4_map_sorted.bam
	mv "$OutputDir"/"$samplename"_IterMap"$iter"/"$samplename"_"$reffile"-iter4_map_sorted.bam.bai "$OutputDir"/"$FinalResults"/"$samplename"_"$reffile"-iter4_map_sorted.bam.bai
	mv "$OutputDir"/"$samplename"_IterMap"$iter"/"$samplename"_IterMap"$iter".log "$OutputDir"/"$FinalResults"/"$samplename"_IterMap"$iter".log
	mkdir "$OutputDir"/"$FinalResults"/BLAST_Hits
	mv "$OutputDir"/*_crunch.txt "$OutputDir"/"$FinalResults"/BLAST_Hits/
	mv "$OutputDir"/"$samplename"_NotHostReads.fasta "$OutputDir"/"$FinalResults"/"$samplename"_NotHostReads.fasta
	mv "$OutputDir"/DCS_Outputs/ "$OutputDir"/"$FinalResults"/DCS_Outputs/
	mv "$OutputDir"/"$FinalResults" "$DIR"/
	
	#Clean up intermediate files
	rm -r "$OutputDir"
	

End=$(date +%s)
TimeTaken=$((End-Start))
echo  | awk -v D=$TimeTaken '{printf "Performed denovoassembly_ONT Analysis in: %02d'h':%02d'm':%02d's'\n",D/(60*60),D%(60*60)/60,D%60}'
