#!/bin/bash
# mmlong assembly
# By Rasmus Kirkegaard
#

# Data:
NANOREADS_FQ=data/nanopore_above4k_filtered.fastq;
ILLUMINAREADS=data/illumina.fastq;

# Settings
THREADS=24;

# PATHS
RACON_PATH=/space/users/rkirke08/Desktop/software/racon/racon/racon/bin/;
SPADES_PATH=/space/users/rkirke08/Desktop/spades/SPAdes-3.10.1-Linux/bin/;
UNICYCLER_PATH=/space/users/rkirke08/Desktop/software/Unicycler/;
MINIASM_PATH=/space/users/rkirke08/Desktop/Miniasm/miniasm/;
MINIMAP_PATH=/space/users/rkirke08/Desktop/Miniasm/minimap/;
MINIMAP2_PATH=/space/users/rkirke08/Desktop/software/minimap2/;

# Prepare directories
mkdir -p temp
mkdir -p results

#####################################################


###################
# Genome assembly #
###################
# Assemble using miniasm (https://github.com/lh3/miniasm)
$MINIMAP_PATH/minimap -Sw5 -L100 -m0 -t $THREADS $NANOREADS_FQ $NANOREADS_FQ | gzip -1 > temp/reads.paf.gz
$MINIASM_PATH/miniasm -f $NANOREADS_FQ temp/reads.paf.gz > temp/miniasm_assembly.gfa
awk '/^S/{print ">"$2"\n"$3}' temp/miniasm_assembly.gfa > temp/miniasm_assembly.fa
# Super fast and produces some promising long contigs

# Polish genome assembly with 2x Racon (https://github.com/isovic/racon)
# Racon1x
$MINIMAP2_PATH/minimap2 -t $THREADS -x map-ont temp/miniasm_assembly.fa $NANOREADS_FQ > temp/mappings1.paf  
$RACON_PATH/racon -t $THREADS $NANOREADS_FQ temp/mappings1.paf temp/miniasm_assembly.fa temp/racon1x.fasta  
# Racon2x
$MINIMAP2_PATH/minimap2 -t $THREADS -x map-ont temp/racon1x.fasta $NANOREADS_FQ > temp/mappings2.paf   
$RACON_PATH/racon -t $THREADS $NANOREADS_FQ temp/mappings2.paf temp/racon1x.fasta temp/racon2x.fasta 
# Clean fasta headers
awk '/^>/{print ">" ++i; next}{print}' temp/racon2x.fasta > results/racon2x_polished.fasta
