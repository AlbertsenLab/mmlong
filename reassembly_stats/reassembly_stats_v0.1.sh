#!/bin/bash

#####################
#### Description ####
#####################

USAGE="$(basename "$0") [-h] [-g fasta -c value -t flag] -- reassembly_stats: Calculation of basic assembly stats.

where:
    -h  Show this help text
    -g  Genome in fasta format
    -c  Number of CPUs to use
    -t  Flag to run pipeline on test data

Run pipeline on test data:
./reassembly_stats_v1.0.sh -t
./reassembly_stats_v1.0.sh -g testdata/genome.fa

Run pipleline on real data:
./reassembly_stats_v1.0.sh -g genome.fa

Requirements - if your programs are named differently, then change the reference in the respective bash script:
quast (version 2.3; http://quast.bioinf.spbau.ru/) 
"
###################
#### Arguments ####
###################

NCPU=2

while getopts ':hz1:g:c:t' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    g) GENOME=$OPTARG;;
    c) NCPU=$OPTARG;;
    t) GENOME=testdata/genome.fa;; 
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

if [ -z ${GENOME+x} ]; then echo "\nERROR: You need to supply a genome e.g. -g genome.fa. See the help information:\n"; echo "$USAGE"; exit 1; fi;

##################
#### Analysis ####
##################

echo "\nCalculating assembly metrics"

### QUAST ###
#quast $GENOME -t 0 > quast_log.txt
#tail -n +6 quast_results/latest/report.txt | awk -F " " '{print $NF}' > quast.txt
#rm -r quast_results

### CheckM ###
#CPATH="$(dirname "$GENOME")"
#checkm lineage_wf -x fa -t $NCPU $CPATH checkm_results -f checkm.txt
cat checkm_results/storage/bin_stats.tree.tsv | tr ',', '\n' | tr ':', '\t' | sed 's/ .//' | sed 's/.\t/\t/' | tail -n +2 > checkm_basestats.txt

cat checkm.txt | sed -n '1~2!p' | sed 's/ \#//g' | sed 's/  /Ø/g' | sed 's/Ø*Ø/\t/g'

### Prokka ###
prokka $GENOME -outdir prokka

#tr '\n', '\t' > res.txt

