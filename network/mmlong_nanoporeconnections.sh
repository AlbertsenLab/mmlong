#!/bin/bash

#####################
#### Description ####
#####################

USAGE="$(basename "$0") [-h] [-threads int -SAM string] -- nanoporeconnections: Produce network file from nanopore sam file.

Requirements:
samtools (v.)
parseTSVtonetwork.R ()
R ()
- dplyr (v. )
"
###################
#### Arguments ####
###################

SAM=/space/users/malber06/2017_long_read/temp/porechopped_filtlong_reads.fastq.gz.minimap.sam;
THREADS=60;
QUALITY=10;

##################
#### Analysis ####
##################

mkdir -p results
mkdir -p temp

# Extract only reads from the sam file that are not flagged as secondary (0x100) or unmapped (0x4) --> (0x104)
echo -e 'READ\tFLAG\tCONTIG' > parsed.tsv
samtools view --threads $THREADS -q $QUALITY -F 0x104  $SAM | cut -d $'\t' -f 1-3 >> temp/parsed.tsv

# Format the connections file using R
Rscript parseTSVtonetwork.R temp/parsed.tsv results/nanopore_network.tsv;
