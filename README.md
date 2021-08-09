# WGS
Scripts for processing whole-genome sequencing data

RefGuidedAlignment
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
  
denovoAssembly
This is a script to perform a de novo assembly of paired-end Illumina NGS data using a given reference database. This is achieved using the following steps: Reads mapping to a given host genome are removed. (The non-host reads are then sub-sampled dependent on number). The non-host reads are the used for a de novo assembly using SPAdes.The contigs output by SPAdes are then BLASTed against a local database (currently databases for influenza virus and avian orthoavulavirus are provided). The top BLAST hit(s) are then used as a reference perform a reference-guided alignment using the raw Illumina reads to produce an intermediate consensus sequence. The intermediate consensus sequence is then used iteratively (iterations = 4) to improve the consensus sequence. A final consensus sequence is then generated using genconsensus.py which utilises the tool iVar (https://github.com/andersen-lab/ivar).

Usage: bash denovoAssembly.sh [options] <Read1.fastq.gz> <Read2.fastq.gz> 
	-G <path to HOST genome File [optional]> 
	-r <path to Reference Database Directory [optional]> 
	-e <Blast e value [optional] default = 0.0001> 
	-t <number of threads/cores to use [optional] default = all>
