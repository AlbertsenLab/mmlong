#!/bin/bash
# mmlong generate coverage
# By Rasmus Kirkegaard
#

#####################
#### Description ####
#####################

USAGE="$(basename "$0") [-h] [-assembly fasta -threads value] -- mmlong: Generate coverage.
where:
    -h  Show this help text
	-assembly Path to assembly.fasta
	-threads Number of threads to user
	illuminasamples file with name(s) of illumina data files
	nanoporesamples file with name(s) of nanopore data files
Requirements - if your programs are named differently, then change the reference in the respective bash script:
minimap2(v. 2.0rc1-r232 https://github.com/lh3/minimap2)
samtools (v. 1.3.1)
bwa (v. <NA> http://bio-bwa.sourceforge.net/)
calc.coverage.in.bam.depth.pl (https://github.com/MadsAlbertsen/multi-metagenome/tree/master/misc.scripts)
"

###################
#### Arguments ####
###################


# Data:
# nanoporesamples file
# illuminasamples file
ASSEMBLY=results/assembly.fasta;

# Settings
THREADS=24;

# PATHS
MINIMAP2_PATH=/space/users/rkirke08/Desktop/software/minimap2/;

# Prepare directories
mkdir -p temp
mkdir -p results

##################
#### Analysis ####
##################

#######################################
# Generate coverage files for binning #
#######################################


echo "Starting to generate nanopore coverage data"
echo "Map reads to assembly using minimap2 (https://github.com/lh3/minimap2)"
while read nanoporesamples
do
NAME=$nanoporesamples;
echo $NAME

# Map nanopore reads to assembly and generate coverage file
$MINIMAP2_PATH/minimap2 -ax map-ont -t $THREADS $ASSEMBLY data/$NAME > temp/$NAME.minimap.sam
# (MA solution from multi metagenome https://github.com/MadsAlbertsen/multi-metagenome/blob/master/misc.scripts/calc.coverage.in.bam.depth.pl)
samtools view --threads $THREADS -Sb  temp/$NAME.minimap.sam -o temp/$NAME.minimap.bam
samtools sort --threads $THREADS temp/$NAME.minimap.bam -o temp/$NAME.minimap.sorted.bam
samtools depth temp/$NAME.minimap.sorted.bam > temp/$NAME.minimap.depth.txt
perl calc.coverage.in.bam.depth.pl -i temp/$NAME.minimap.depth.txt -o results/$NAME.minimap.coverage.csv

done < nanoporesamples


echo "Starting to generate illumina coverage data"
echo "Map reads to assembly using bwa (http://bio-bwa.sourceforge.net/)"
bwa index $ASSEMBLY
while read illuminasamples
do
NAME=$illuminasamples;
echo $NAME
bwa mem -t $THREADS $ASSEMBLY data/$NAME > temp/$NAME.aln.sam

echo "Generate coverage data"
# (MA solution from multi metagenome https://github.com/MadsAlbertsen/multi-metagenome/blob/master/misc.scripts/calc.coverage.in.bam.depth.pl)
samtools view --threads $THREADS -Sb  temp/$NAME.aln.sam  >  temp/$NAME.aln.bam
samtools sort  --threads $THREADS temp/$NAME.aln.bam -o temp/$NAME.aln.sorted.bam
samtools depth temp/$NAME.aln.sorted.bam > temp/$NAME.depth.txt
perl calc.coverage.in.bam.depth.pl -i temp/$NAME.depth.txt -o results/$NAME.bwa.coverage.csv
done < illuminasamples


#############################################
# Perform binning of the assembly elsewhere #
# Using the mmgenome package or             #
# an automated binning solution             #
#############################################