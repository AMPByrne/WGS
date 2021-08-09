denovoAssembly

Overview
This is a script to perform a de novo assembly of paired-end Illumina NGS data using a given reference database.

This is achieved using the following steps: Reads mapping to a given host genome are removed. (The non-host reads are then sub-sampled dependent on number).

The non-host reads are the used for a de novo assembly using SPAdes.

The contigs output by SPAdes are then BLASTed against a local database (currently databases for influenza virus and avian orthoavulavirus are provided).

The top BLAST hit(s) are then used as a reference perform a reference-guided alignment using the raw Illumina reads to produce an intermediate consensus sequence.

The intermediate consensus sequence is then used iteratively (iterations = 4) to improve the consensus sequence.

A final consensus sequence is then generated using genconsensus.py which utilises the tool iVar (https://github.com/andersen-lab/ivar).

Version History
version 0 - The initial version of this script was based on FluSeqID (https://github.com/ellisrichardj/FluSeqID) but was re-written by James Seekings to run on CentOS rather than Biolinux.

verion 1.0 - Velvet was replaced with SPAdes due to ongoing issues. This script was then renamed to denovoAssembly.sh as it had changed significantly form the original FluSeqID script.

version 2.0 - Removed the use vcf2consensus.pl for generating the final consensus. Instead this has now been replaced with genconsensus.py, which uses the program iVar to generate the final consensus. Dependencies updated accordingly. Fixed the inclusion of the data in the results directory.

version 3.0 - Updated to work with genconsensus.py version 2.0, which allows the user to define variables for iVar. The usage function has been updated accordingly. The flag for threads has been changed from -t to -T to avoid conflict with iVar flags. Note - in version 2.0 of genconsensus.py, the default minimum frequency threshold has been changed from 0.8 to 0, to ensure the base that is observed in the majority of reads (accounting for read quality) is used to call the consensus. The default minimum read depth is still 1, and the default character used for nucleotides without coverage is still -.

Dependencies
BWA samtools bamtools SPAdes blast+ picard GATK iVar bcftools genconsensus.py vcf2consensus.pl

Installation
This repo contains the denovoAssembly-v2.yml file to create an env in conda with all the appropriate software tools installed. To create this env use the following commands:

cp /home/$USER/mnt/VI6Bioinformatics/Central_Pipelines/denovoAssembly/denovoAssembly-v2.yml /home/$USER/

conda env create -f /home/$USER/denovoAssembly-v2.yml
The conda env can then be activated using the following command:

conda activate denovoAssembly-v2

NOTE: you need to activate this env before trying to run the script.

The conda env contains all the dependencies, except: vcf2consensus.pl, genconsensus.py and gatk.py scripts, which are included in /home/$USER/mnt/VI6Bioinformatics/Cental_Pipelines/. However, when using this script on the VI6 SCEv3 architecture you won't need to install or move these.

On the new SCEv3 virtual mahcines, you need to install java. You can do this by following these commands:

sudo apt update 
NOTE: This will tell you if any packages can be ugraded - can be skipped

sudo apt upgrade 
NOTE: This will install updates/upgrades

sudo apt install default-jdk
There are also issues currently with the version of openssl that samtools installs which causes issues with bcftools. To get around this, you'll have to install bcftools directly on your machine until this is addressed. To do this run the following command:

sudo apt-get install bcftools
Usage
Once the conda env is activated the script is used with the following command changing the inputs (NOTE: the path to the script is quite long, so you can copy and paste it into a terminal):

bash /home/$USER/mnt/VI6Bioinformatics/Central_Pipelines/denovoAssembly/denovoAssembly.sh [options] <Read1.fastq.gz> <Read2.fastq.gz>
-G <path to HOST genome File  [optional] default = Gallus gallus genome>
-r <path to Reference Database Directory [optional] default = Influenza database>
-e <Blast e value [optional] default = 0.0001>
-t <number of threads/cores to use [optional] default = all>
NOTE: the path to the script is quite long, so you may want to make a symbolic link to it - see the tutorial for making symbolic links in /home/alexbyrne/mnt/VI6Bioinformatics/Central_Pipelines/General_Tutorials/

Available Viral Reference Databases:

Influenza - default

AOAV1 - use -r home/$USER/mnt/VI6Bioinformatics/Central_Pipelines/denovoAssembly/Viral_Reference_Databases/AOAV_Database

Hepeviridae - use -r /home/$USER/mnt/VI6Bioinformatics/Central_Pipelines/denovoAssembly/Viral_Reference_Databases/Hepeviridae_Database

Paramyxoviridae - use -r /home/$USER/mnt/VI6Bioinformatics/Central_Pipelines/denovoAssembly/Viral_Reference_Databases/Paramyxoviridae_Database

Canine distemper viru - use -r /home/$USER/mnt/VI6Bioinformatics/Central_Pipelines/denovoAssembly/Viral_Reference_Databases/CanineDistemper_Database
NOTE: the Influenza Database is kept offline and is for INTERNAL USE ONLY and NOT FOR DISTRIBUTION. You need a GISAID account (https://www.gisaid.org/) to be able to use this database to comply with GISAID's Terms and Conditions. Do not distribute this database to any person(s) who does not have a GISAID account!

Available Host Genomes: Chicken - default

Human - use -G /home/$USER/mnt/VI6Bioinformatics/Central_Pipelines/denovoAssembly/Host_Genomes/Human_Genome/hg38.fa.gz

Ferret - use -G /home/$USER/mnt/VI6Bioinformatics/Central_Pipelines/denovoAssembly/Host_Genomes/MustelaPutoriusFuroGenome/AEYP01.fasta.gz

Pig - use -G /home/$USER/mnt/VI6Bioinformatics/Central_Pipelines/denovoAssembly/Host_Genomes/SusScrofaGenome/GCF_000003025.6_Sscrofa11.1_genomic.fna.gz

Fox  -use -G /home/$USER/mnt/VI6Bioinformatics/Central_Pipelines/denovoAssembly/Host_Genomes/VulpesVulpesGenome/GCF_003160815.1_VulVul2.2_genomic.fna.gz

Gray seal -se -G /home/$USER/mnt/VI6Bioinformatics/Central_Pipelines/denovoAssembly/Host_Genomes/HalichoerusGrypusGenome/GCA_012393455.1_Tufts_HGry_1.1_genomic.fna.gz
The host genomes are too large to be hosted on Github, so are only kept locally.

Acknowledgements

James Seekings (APHA) Alex Byrne (APHA) Ben Mollett (APHA)

