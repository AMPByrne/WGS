#!/bin/bash
set -e

#Reference-guided alignment of paired-end Illumina reads.
#version 1.0 - initial script to map paired-end Illumina reads to a reference sequence.
#version 2.0 - added removal of PCR duplicates using Samtools fixmate and rmdup, also added a function to generate mapping stats using Samtools flagstat pre- and post-duplicate removal so you can monitor the number of PCR duplicates and how the removal of them affects your mapping quality. Also added a timer to show how long the entire process took. Added requisite programs/tools as well.
#
#input: RefGuidedAlignment.sh Samplename reference.fasta R1reads.fastq.gz R2reads.fastq.gz outputlocation numberofcores

Samplename="$1"
reference="$2"
R1reads="$3"
R2reads="$4"
Outputlocation="$5"
Cores="$6"

if [ $# -lt 5 ]; then
echo "
input: RefGuidedAlignment.sh Samplename reference.fasta R1reads.fastq.gz R2reads.fastq.gz outputlocation numberofcores
"
exit 1
fi
#
#Required programs/tools:
#BWA
#Samtools
#bcftools
#vcf2consensus.pl
#time the process
Start=$(date +%s)
#
#make a new results directory and cd to this directory. replace [directory] with the directory you wish to use.
mkdir /[directory]/"$Samplename"_Temp_Results
cd /[directory]/"$Samplename"_Temp_Results
#
#Index reference genome
echo Indexing reference sequence
bwa index "$reference"
#
#align paired end reads to the indexed reference and output a bam file
echo Aligning reads to reference sequence
bwa mem -t"$Cores" "$reference" "$R1reads" "$R2reads" > "$Samplename".sam
#
#convert the mapped sam file to a bam file and fill in mate coordinates, ISIZE and mate related flags
echo Converting sam to bam and filling in mate coordinates
samtools fixmate -O bam "$Samplename".sam "$Samplename".bam
#
#convert the mapped sam file to a bam file
#echo Converting sam to bam
#samtools view -b -@"$Cores" "$Samplename".sam -o #"$Samplename".bam
#
#generate mapping stats
echo Generating mapping stats
samtools flagstat "$Samplename".bam > "$Samplename".bam.STATS.pre.dup.removal.txt
#
#remove duplicate reads
echo Removing duplicate reads
samtools rmdup -sS "$Samplename".bam "$Samplename".rmdup.bam
#
#generate mapping stats post duplicate removal
echo Generating mapping stats post duplicate removal
samtools flagstat "$Samplename".rmdup.bam > "$Samplename".bam.STATS.rmdup.txt
#
#sort the bam file
echo Sorting bam file
samtools sort "$Samplename".rmdup.bam -o "$Samplename".rmdup.sorted.bam
#
#index the bam file
echo Indexing bam file
samtools index -b "$Samplename".rmdup.sorted.bam
#
#create a pileupfile of the bwa mapping
echo Creating mpileup file
samtools mpileup -L 10000 -Q 1 -AEupf "$reference" "$Samplename".rmdup.sorted.bam >  "$Samplename".rmdup.sorted.bam.mpileup
#
#call the varients between the data and the reference and create a vcf file
echo Calling varients
bcftools call -c "$Samplename".rmdup.sorted.bam.mpileup > "$Samplename".rmdup.sorted.bam.vcf
#
#create a new consensus sequence
echo Generating new consensus sequence
vcf2consensus.pl consensus -f "$reference" "$Samplename".rmdup.sorted.bam.vcf > "$Samplename".fasta
#
#[optional] you then then do some stats on the mapping
#echo Generating mapping stats
#samtools flagstat "$Samplename".bam > "$Samplename"_STATS.txt
#
#just clearing up now and move the final files to the output location
echo Time to clean up
rm "$Samplename".sam
rm "$Samplename".bam
rm "$Samplename".rmdup.bam
rm "$Samplename".rmdup.sorted.bam.mpileup
rm "$Samplename".rmdup.sorted.bam.vcf
cd "$Outputlocation"
mkdir "$Samplename"_Results
cp "$reference".* "$Outputlocation"/"$Samplename"_Results
rm "$reference".*
mkdir "$Outputlocation"/"$Samplename"_Results/Reads
mv "$R1reads" "$Outputlocation"/"$Samplename"_Results/Reads
mv "$R2reads" "$Outputlocation"/"$Samplename"_Results/Reads
#replace [directory] with the directory you wish to use.
mv /[directory]/"$Samplename"_Temp_Results/* "$Outputlocation"/"$Samplename"_Results
cd /[directory]/
rm -r "$Samplename"_Temp_Results
End=$(date +%s)
TimeTaken=$((End-Start))
echo "All finished,you're results are here:' $Outputlocation"
echo | awk -v D=$TimeTaken '{printf "RefGuidedAlignment.sh took %02d'h':%02d'm':%02d's'\n",D/(60*60),D%(60*60)/60,D%60}'
