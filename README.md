# WGS
Scripts for processing whole-genome sequencing data

## RefGuidedAlignment
The script RefGuidedAlignment_Public.sh uses BWA/Minimap2 samtools and iVar to map raw NGS reads (from either Illumina or Oxford Nanopore sequencers) to a reference sequence in a fasta file.

Usage: bash RefGuidedAlignment_Public.sh [options] <read(s).fastq>
	-r path to reference genome file [Required]
	-s sequencer used, Illumina or ONT [Required]
	-o path to output location [optional] default = location of read file(s)
	-c cores [optional] default = all
	-t the minimum frequency at which a nucleotide will be used to call the consensus. INT 0 - 1.[optional] deault = $MinFreq
	-m minimum read depth at which iVar call a consensus [optional] default = 1
	-n (N or -) character to print if less than minimum coverage (-m) [optional] default = -
  
genconsensus.py is a dependency of RefGuidedAlignment_Public.sh and denovoAssembly_Public.sh and converts the sorted bam file and reference into a final consensus sequence.
  
The RefGuidedAlignment-v2.yml file is a conda env file containing all the software/tools required to run RefGuidedAlignment.sh and genconsensus.py.
  
## denovoAssembly
There are two scripts to perform a hybrid de novo assembly/reference-guided alignment of raw Illumina paired-end (denovoAssembly_Public.sh) or Oxford Nanopore reads (denovoAssembly_ONT_Public.sh) using a given the appropriate host genome and a reference database. 

This is achieved using the following steps: 
1. Reads mapping to a given host genome are removed. The non-host reads are then sub-sampled dependent on number. For ONT data, the reads are also aligned against the DNA Control Sample sequence (DNA Lambda phage) that is provided with the V14 kits. This is useful for troubleshooting.
2. The non-host reads are the used for a de novo assembly using SPAdes. (Illumina data only).
3. The contigs output by SPAdes (Illumina data only) or the reads themselves (ONT data) are then BLASTed against a local database.
4. The top BLAST hit(s) are then used as a reference perform a reference-guided alignment using the raw reads to produce an intermediate consensus sequence. The intermediate consensus sequence is then used iteratively (iterations = 4) to improve the consensus sequence.
5. A final consensus sequence is then generated. For Illumina data this uses genconsensus.py which utilises the tool iVar (https://github.com/andersen-lab/ivar). For ONT data, medaka is used to generate the intermediate and final consensus, with the later masked based on the read depth (min read of 10).


Usage: bash denovoAssembly_Public.sh [options] <Read1.fastq.gz> <Read2.fastq.gz> 
	-G <path to HOST genome File> 
	-r <path to Reference Database Directory> 
	-e <Blast e value [optional] default = 0.0001> 
	-t <number of threads/cores to use [optional] default = all>

 Usage: bash denovoAssembly_ONT_Public.sh  [options] -G <HostGenome> -r <ReferenceDatabase> <Reads.fastq.gz> 
	help Displays this usage
	-G <path to HOST genome File>
	-r <path to Reference Database Directory>
	-e <Blast e value [optional] default = 0.0001>
	-T <number of threads/cores to use [optional] default = all>
	-m <minimum read depth at which a consensus sequence is called. Region with read coverage below this will be masked [optional] default = 10>

The denovoAssembly-v2.yml file is a conda env file containing all the software/tools required to run denovoAssembly_Public.sh and genconsensus.py.

The denovoAssembly-v3.yml file is a conda env file containing all the software/tools required to run denovoAssembly_ONT_Public.sh.
