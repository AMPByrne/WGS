#!/bin/bash
set -e

#Usage: find_gaps_degenerates.sh -c consensus.fasta -b bamfile

while getopts 'c:b:' opt ; do
    case $opt in
        c) PathToConsensusFile=$OPTARG ;;
        b) PathToBamFile=$OPTARG ;;
    esac
done

if [ $# -lt 2 ]; then
    echo "Usage: find_gaps_degenerates.sh -c consensus.fasta -b bamfile"
fi

pysamstats --type variation --fasta "$PathToConsensusFile" "$PathToBamFile" > "$PathToBamFile".mismatch.txt

python bin/filter_pysamstatsfile.py "$PathToBamFile".mismatch.txt
