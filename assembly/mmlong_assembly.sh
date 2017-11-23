#!/bin/bash
# mmlong assembly
# By Rasmus Kirkegaard
#
#####################
#### Description ####
#####################

USAGE="$(basename "$0") [-h] [-nanopore fastq -threads value -minlength value] -- mmlong: Run long read assembly.
where:
    -h  Show this help text
	-nanopore Path to your fastq file
	-threads Number of threads to user
	-minlength Minimum read length requirement
Requirements - if your programs are named differently, then change the reference in the respective bash script:
minimap (v. 0.2-r124-dirty https://github.com/lh3/minimap)
miniasm (v. 0.2-r168-dirty https://github.com/lh3/miniasm)
minimap2(v. 2.0rc1-r232 https://github.com/lh3/minimap2)
racon (v. <NA> https://github.com/isovic/racon)
"

###################
#### Arguments ####
###################

# Data:
NANOREADS_FQ=data/nanopore_above4k_filtered.fastq;

# Settings
THREADS=24;
min_length=4000;

# PATHS
RACON_PATH=/space/users/rkirke08/Desktop/software/racon/racon/racon/bin/;
MINIASM_PATH=/space/users/rkirke08/Desktop/Miniasm/miniasm/;
MINIMAP_PATH=/space/users/rkirke08/Desktop/Miniasm/minimap/;
MINIMAP2_PATH=/space/users/rkirke08/Desktop/software/minimap2/;

##################
#### Analysis ####
##################

# Prepare directories
mkdir -p temp
mkdir -p results

#############################
# Trim adapters w. Porechop #
#############################
source porechop.env.sh add
porechop -i $NANOREADS_FQ -o temp/porechopped_reads.fastq --threads $THREADS --min_split_read_size $min_length
source porechop.env.sh remove

########################################
# Filter low quality reads w. Filtlong #
########################################

# Without reference
filtlong --min_length $min_length --keep_percent 90 temp/porechopped_reads.fastq | gzip > temp/porechopped_filtlong_reads.fastq.gz

NANOREADS_FQ_HQ=temp/porechopped_filtlong_reads.fastq.gz;

###################
# Genome assembly #
###################
# Assemble using miniasm (https://github.com/lh3/miniasm)
$MINIMAP_PATH/minimap -Sw5 -L100 -m0 -t $THREADS $NANOREADS_FQ_HQ $NANOREADS_FQ_HQ | gzip -1 > temp/reads.paf.gz
$MINIASM_PATH/miniasm -f $NANOREADS_FQ_HQ temp/reads.paf.gz > temp/miniasm_assembly.gfa
awk '/^S/{print ">"$2"\n"$3}' temp/miniasm_assembly.gfa > temp/miniasm_assembly.fa
# Super fast and produces some promising long contigs

# Polish genome assembly with 2x Racon (https://github.com/isovic/racon)
# Racon1x
$MINIMAP2_PATH/minimap2 -t $THREADS -x map-ont temp/miniasm_assembly.fa $NANOREADS_FQ_HQ > temp/mappings1.paf  
$RACON_PATH/racon -t $THREADS $NANOREADS_FQ_HQ temp/mappings1.paf temp/miniasm_assembly.fa temp/racon1x.fasta  
# Racon2x
$MINIMAP2_PATH/minimap2 -t $THREADS -x map-ont temp/racon1x.fasta $NANOREADS_FQ_HQ > temp/mappings2.paf   
$RACON_PATH/racon -t $THREADS $NANOREADS_FQ_HQ temp/mappings2.paf temp/racon1x.fasta temp/racon2x.fasta 
# Clean fasta headers
awk '/^>/{print ">" ++i; next}{print}' temp/racon2x.fasta > results/assembly.fasta
