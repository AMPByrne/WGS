'''
Alex Byrne Jul 2020
Animal and Plant Health Agency

This is a script to generate a consensus sequence:
1. An mpileup file is generated from a reference fasta file and an indexed, sorted bam file using samtools. 
2. The mpileup file is then used to generate a consensus sequence using iVar. 

iVar is a software tool that allows you to specify the proportions at which a degenerate nucleotide is called, and to also set a minimum depth for which a consensus will be called. If the depth is below this threshold it will insert a defined character (N or -) at these positions. For more information visit the iVar manual: https://andersen-lab.github.io/ivar/html/manualpage.html 
The only problem is that iVar currently only works with reference fasta files contianing a single sequence. This would be sufficient for non-segmented references (bacterial chromosomes, non-segmented viruses), however, since this will primarily be used for influenza sequences, I have had to come up with this script as a work around. After speaking to the creator of iVar, they will implement the ability to process reference fasta files containing multiple sequences, this may take some time. 
This script will work for segmented/non-segmented reference fasta files, however, if you are not using a segmented reference, the following command will produce the same output:
        samtools mpileup -aa -A -d 99999 -Q 0 -B -f /path/to/reference.fasta, path/to/bam | ivar consensus -t minfrequencythreshold -m mindepth -n- -p consensus_prefix

The default options for iVar are:
        -t 0 - this will call the majority or most common base as the consensus.See https://andersen-lab.github.io/ivar/html/manualpage.html#autotoc_md19 for more info.
        -m 1 - minimum depth to call a consensus is 1 read.
        -n- - any positions with depth < the integer provided for -m will be replace with -

Version 2.0 - enabled the ability for the user to change the default options for iVar. Renamed the variables to be more similar to how they are described in the iVar manual. Changed the default minimum frequency threshold from 0.7 to 0 so that the majority or most common base is called as a consensus - will reduce the number degenerates called.

Dependencies:
    Samtools
    iVar
      
'''
import os
import subprocess
import sys
import getopt
'''
reference = '/home/p992561/mnt/Alex/Consensus_Testing/ND-Test_Results/AJ880277.1 Pigeon paramyxovirus-1.fas'
bam = '/home/p992561/mnt/Alex/Consensus_Testing/ND-Test_Results/ND-Test.rmdup.sorted.bam'
'''
VERSION = 2.0
#default iVar options
MinFreq = '0'
MinDepth = '1'
NoCoverage = '-'


def usage():
    print("genconsensus.py VERSION:",VERSION,
         "\n\tUsage: python",sys.argv[0],"[Options] -r <reference.fasta> -b <sortedbam.bam>"
          "\n\t-r, --reference reference fasta file [Required]"
          "\n\t-b, --bam sorted bam file [Requried]"
          "\n\t-t, --minfreq INT [0-1] the minimum frequency at which a nucleotide will be used to call the consensus. [OPTIONAL] DEFAULT =",MinFreq, 
          "\n\t-m, --mindepth INT the minimum read depth at which a consensus base is called. [OPTIONAL] DEFAULT =",MinDepth,
          "\n\t-n, --nocoverage -/N any positions with read depth less than -m will be given the - or N character. iVar only allows the use of one of these two characters. [OPTIONAL] DEFAULT =",NoCoverage,
          "\n\t-h, --help") 

#get user options
try:
    opts, args = getopt.getopt(sys.argv[1:], 'r:b:t:m:n:h', ["reference=", "bam=", "minfreq=", "mindepth=", "nocoverage=", "help" ])
except getopt.GetoptError:
    usage()
    sys.exit(1)
    
for opt, arg in opts:
    if opt in ("-h", "--help"):
        usage()
        sys.exit(1)
    elif opt in ("-r", "--reference"):
        reference = arg
    elif opt in ("-b", "--bam"):
        bam = arg
    elif opt in ("-t", "--minfreq"):
        MinFreq = arg
    elif opt in ("-m", "--mindepth"):
        MinDepth = arg
    elif opt in ("-n", "--nocoverage"):
        NoCoverage = arg
        
if not opts:
    print("*****Error - A reference fasta file and a saorted bam file must be provided*****")
    usage()
    sys.exit(1)

    
#index the reference to allow you to generate an mpileup for specific regions within the reference
print("Indexing the reference sequence using samtools faidx")
subprocess.run(["samtools", "faidx", reference])

#Find the fasta headers within the reference and assign the variable header

for line in open(reference).readlines():
    if line.startswith(">"):
        header = line
        #To generate a region for samtools to make the mpileup file, remove the > and assign to variable region
        region = header[1:].strip()
        #To generate a filename we can use for the output files, remove any special characters commonly found in fasta headers (/ and |)
        outputname = region.replace("/","~~")
        outputname = outputname.replace("|","~")
        print("Generating consensus for region: ", region) 
        #Generate mpilup for each region and pass to iVar to make the consensus sequence for each region
        samtools = subprocess.run(["samtools", "mpileup", "-aa", "-A", "-d 99999", "-Q 0", "-B", "-r", region, "-f", reference, bam], stdout=subprocess.PIPE)
        subprocess.run(["ivar", "consensus", "-t", MinFreq, "-m", MinDepth, "-n", NoCoverage, "-p", outputname], input=samtools.stdout)

#Now we're going to merge the consensus sequences that were generated by iVar.
print("Merging the separate consensus files")
#Create an output file called final_consensus.fasta       
with open("final_consensus.fasta", "w") as outfile:
    #Find all the .fa files generate by iVar and merge into final_consensus.fasta
    for file in os.listdir():
        if file.endswith(".fa"):
            with open(file, "r") as readfile:
                outfile.write(readfile.read() + "\n\n")

#Just going to reformat the headers now back to how they were in the reference
with open("final_consensus.fasta", "r") as final_consensus_fasta:
    lines = final_consensus_fasta.readlines()
with open("final_consensus.fasta", "w") as f:
    for line in lines:
        if not line.isspace():
            line = line.replace("~~","/")
            line = line.replace("~","|")
            line = line.replace(">Consensus_", ">")
            f.write(line)

#Remove the qual.txt files and the individual consensus files
directory = os.listdir()
for file in directory:
    if file.endswith(".qual.txt"):
        os.remove(file)
    #if file.endswith(".fa"):
     #   os.remove(file)
        
print("Consensus sequence generated")
       
    

                
        
        
        
     
