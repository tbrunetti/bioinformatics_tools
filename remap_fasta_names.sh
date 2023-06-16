#!/bin/bash

#BSUB -q normal # name of queue to use
#BSUB -J mouseGRCm38.p6_multiStrain # job name
#BSUB -P general
#BSUB -W 15:00 # max walltime of 15 hours
#BSUB -e mouseGRCm38.p6_multiStrain_02102021.err # error file
#BSUB -o mouseGRCm38.p6_multiStrain_02102021.log # stdout file
#BSUB -n 10 # number of CPU cores
#BSUB -R "span[hosts=1]" # span the 5 CPU cores requested above across 1 node
#BSUB -R "rusage[mem=5GB]" # 5BG per core/slot is required
#BSUB -M 5G # kill the job if 5GB is exceeded per core/slot

while read line; do
	#echo "$line"
	findValue="$(cut -d' ' -f1 <<< $line)"
	replaceValue="$(cut -d' ' -f2 <<< $line)"
	sed -i "s?$findValue?$replaceValue?" GCF_000001635.26_GRCm38.p6_genomic_chromNamesRemapped_OtherStrainsIncluded.fasta
done <GRCm38.p6_chromNamesRemapped_OtherStrainsIncluded.txt

