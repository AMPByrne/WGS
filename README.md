# WGS
Scripts for processing whole-genome sequencing data

The script RefGuidedAlignment uses BWA, samtools and iVar to map raw NGS reads (from either Illumina or Oxford Nanopore sequencers) to a reference sequence in a fasta file.

Usage: bash RefGuidedAlignment.sh [options] <read(s).fastq>
	-r path to reference genome file [Required]
	-s sequencer used, Illumina or ONT [Required]
	-o path to output location [optional] default = location of read file(s)
	-c cores [optional] default = all
	-t the minimum frequency at which a nucleotide will be used to call the consensus. INT 0 - 1.[optional] deault = $MinFreq
	-m minimum read depth at which iVar call a consensus [optional] default = $MinDepth
	-n (N or -) character to print if less than minimum coverage (-m) [optional] default = $NoCoverage
  
  genconsensus.py is a dependency of RefGuidedAlignment.sh and converts the sorted bam file and reference into a final consensus sequence.
  
  The RefGuidedAlignment-v2.yml file is a conda env file containing all the software/tools required to run RefGuidedAlignment.sh and genconsensus.py.
