#!/bin/bash
# mmlong reassembly workflow
# By Rasmus Kirkegaard and Søren Karst
VERSION=1.0.4

################################################################################
### Preparation ----------------------------------------------------------------
################################################################################

### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [-n file -i file -b folder -d folder -m folder -t value -x] 
-- mmlong reassembly workflow v. $VERSION: Based on binned scaffolds recruit reads and perform Unicycler reassembly.

where:
    -h  Show this help text.
    -n  List of nanopore datafiles to be used for re-assembly. Default: ./np_asmb.txt.
    -i  List of Illumina datafiles to be used for re-assembly. Default: ./ilm_asmb.txt.
    -b  Folder containing bins in fasta format. Default: ./binning/bins .
    -d  Folder containing trimmed reads. Default: ./trimmed_data .
    -m  Folder containing mapped reads. Default: ./metagenome_mapping .
    -t  Number of threads to use. Default: 10.
    -x  Run pipeline on test data.

Requirements:
- Test

"
### Customizable Arguments -----------------------------------------------------

# Paths to dependencies
SAMTOOLS=/space/users/smk/Software/SAMtools/bin/samtools
RACON=/space/users/rkirke08/Desktop/software/racon/racon/racon/bin/racon;
SPADES=/space/users/smk/Software/Spades3.11/SPAdes-3.11.1-Linux/bin/spades.py;
UNICYCLER=/space/users/smk/bin/unicycler;
FILTLONG=/space/sharedbin/bin/filtlong;
QUAST=/space/sharedbin/bin/quast;

### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hzn:i:b:d:m:t:x' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    n) NP_ASMB=$OPTARG;;
    i) ILM_ASMB=$OPTARG;;
    b) BIN_DIR=$OPTARG;;
    d) TRIM_DIR=$OPTARG;;
    m) MAP_DIR=$OPTARG;;
    t) THREADS=$OPTARG;;
    x) NP_ASMB=`cat ./np_asmb.txt`;ILM_ASMB=`cat ./ilm_asmb.txt`;\
       BIN_DIR=./binning/bins; TRIM_DIR=./trimmed_data; MAP_DIR=./metagenome_mapping;\
       THREADS=40;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${THREADS+x} ]; then THREADS=10; fi;

### Data names
NP_ASMB_NAME=`echo "$NP_ASMB" | sed 's/\..*$//g'`
ILM_ASMB_NAME=`echo "$ILM_ASMB" | sed 's/[12]\..*$//g' | sort | uniq`

################################################################################
### Workflow -------------------------------------------------------------------
################################################################################

### Read recruitment -----------------------------------------------------------
if [ ! -d "reassembly" ]; then
mkdir reassembly

# Prepare bin folders
for BIN_PATH in $BIN_DIR/*.fa*
do
  BIN_NAME=`basename "$BIN_PATH" | sed 's/\..*$//g'`
  mkdir reassembly/$BIN_NAME
  cp $BIN_PATH reassembly/$BIN_NAME/${BIN_NAME}_massembly.fa
  grep ">" $BIN_PATH | tr -d ">" |  tr -d '\15\32' \
  > reassembly/$BIN_NAME/scaffold.list
  mkdir reassembly/$BIN_NAME/data
done

# Extract reads
for BIN_PATH in reassembly/*
do
  for ILM_NAME in $ILM_ASMB_NAME
  do
    IR=`$SAMTOOLS view -@ $THREADS -q 10 -F 0x104 $MAP_DIR/${ILM_NAME}_cov.bam |\
    awk 'NR==FNR{a[$0];next}$3 in a{print $1}' $BIN_PATH/scaffold.list - |\
    sort | uniq` # Subsetting disabled
    LC_ALL=C grep -Fwf <(echo "$IR") -A 3 $TRIM_DIR/${ILM_NAME}1_trim.fq | \
    sed '/^--$/d' > $BIN_PATH/data/${ILM_NAME}1_trim_subset.fq &
    LC_ALL=C grep -Fwf <(echo "$IR") -A 3 $TRIM_DIR/${ILM_NAME}2_trim.fq | \
    sed '/^--$/d' > $BIN_PATH/data/${ILM_NAME}2_trim_subset.fq    
  done 
done

for BIN_PATH in reassembly/*
do
  for NP_NAME in $NP_ASMB_NAME
  do
    IR=`$SAMTOOLS view -@ $THREADS -q 10 -F 0x104 $MAP_DIR/${NP_NAME}_cov.bam |\
    awk 'NR==FNR{a[$0];next}$3 in a{print $1}' $BIN_PATH/scaffold.list - |\
    sort | uniq`
    LC_ALL=C grep -Fwf <(echo "$IR") -A 3 $TRIM_DIR/${NP_NAME}_trim.fq | \
    sed '/^--$/d' > $BIN_PATH/data/${NP_NAME}_trim_subset.tmp  
    #$FILTLONG --target_bases 150000000 \
    #--length_weight 10 $BIN_PATH/data/${NP_NAME}_trim_subset.tmp \ #Subsetting disabled
    cat $BIN_PATH/data/${NP_NAME}_trim_subset.tmp > $BIN_PATH/data/${NP_NAME}_trim_subset.fq
    rm $BIN_PATH/data/${NP_NAME}_trim_subset.tmp
  done 
done
fi

### Reassembly -----------------------------------------------------------------
for BIN_PATH in reassembly/*
do
  ASMBS="$BIN_PATH/*_massembly.fa"
  for ILM_NAME in $ILM_ASMB_NAME
  do
    for NP_NAME in $NP_ASMB_NAME
    do
      ASMB_NAME=${NP_NAME}+${ILM_NAME}
      ASMBS="$ASMBS $BIN_PATH/${ASMB_NAME}_assembly.fa"
      if [ ! -d "$BIN_PATH/$ASMB_NAME" ]; then
        mkdir $BIN_PATH/$ASMB_NAME
        $UNICYCLER -1 $BIN_PATH/data/${ILM_NAME}1_trim_subset.fq \
        -2 $BIN_PATH/data/${ILM_NAME}2_trim_subset.fq \
        -l $BIN_PATH/data/${NP_NAME}_trim_subset.fq --threads $THREADS \
        --spades_path $SPADES --no_correct --min_kmer_frac 0.3 --kmer_count 3 \
        --no_pilon --racon_path $RACON --keep 3 -o $BIN_PATH/$ASMB_NAME
        cp $BIN_PATH/$ASMB_NAME/assembly.fasta \
        $BIN_PATH/${ASMB_NAME}_assembly.fa
      fi       
    done
  done
  if [ ! -d "$BIN_PATH/quast_stats" ]; then
    $QUAST $ASMBS -o $BIN_PATH/quast_stats -R $BIN_PATH/*_massembly.fa
  fi
done

### Logfile --------------------------------------------------------------------
LOG_NAME="mmlong_reassembly_log_$(date +%s).txt"

echo "Input arguments:
mmlong reassembly workflow script version: $VERSION
Path to genome bins: $BIN_DIR
Nanopore data list for assembly: $NP_ASMB
Illumina data list for assembly: $ILM_ASMB
Path to trimmed data: $TRIM_DIR
Path to mapped data: $MAP_DIR
Number of threads: $THREADS
Nanopore min read length: $NP_MINLENGTH" >> $LOG_NAME

################################################################################
### Exit -----------------------------------------------------------------------
################################################################################
exit 0

################################################################################
### Testing --------------------------------------------------------------------
################################################################################
NP_ASMB=`cat ./np_asmb.txt`;
ILM_ASMB=`cat ./ilm_asmb.txt`;
BIN_DIR=./binning/bins;
TRIM_DIR=./trimmed_data;
MAP_DIR=./metagenome_mapping;
THREADS=40;
NP_ASMB_NAME=`echo "$NP_ASMB" | sed 's/\..*$//g'`
ILM_ASMB_NAME=`echo "$ILM_ASMB" | sed 's/[12]\..*$//g' | sort | uniq`


