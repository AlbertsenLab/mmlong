#!/bin/bash
# mmlong metagenome data generation
# By Rasmus Kirkegaard and Søren Karst
# Version
VERSION=1.1.4

################################################################################
### Preparation ----------------------------------------------------------------
################################################################################

### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [-d path -n file -m file -i file -f -t value -l value -x -k] 
-- mmlong metagenome data generation v. $VERSION: Read filtering, metagenome assembly, 
read coverage estimation, taxonomic classification and detection of SSU rRNA.

where:
    -h  Show this help text.
    -d  Sequencing data folder.
    -n  List of Nanopore data files to be used for assembly.
    -m  List of Nanopore data files to be used for read coverage estimates.
    -i  List of Illumina PE data files to be used for read coverage estimates.
        File names should end with \"1.fq/2.fq\" or \"1.fastq/2.fastq\". 
    -f  Sequencing data is pre-trimmed.
    -t  Number of threads to use.
    -l  Minimum nanopore reads length.
    -x  Run pipeline on test data.
    -k  Keep temporary files.

Run pipeline on test data:
./mmlong_metagenome_workflow_v1.sh -x
./mmlong_metagenome_workflow_v1.sh -d data -n np_asmb.txt -m np_cov.txt
-i ilm_cov.txt -t 40 -l 8000

Requirements:

## Write dependencies and versions ##
"

### Customizable Arguments -----------------------------------------------------

# Adaptors
NEX_ADP1=CTGTCTCTTATACACATCT # Illumina Nextera adaptor sequences
NEX_ADP2=CTGTCTCTTATACACATCT # Illumina Nextera adaptor sequences
TRU_ADP1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA # Illumina TruSeq and NEB Nebnext adaptor sequences
TRU_ADP2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT # Illumina TruSeq and NEB Nebnext adaptor sequences

# Paths to dependencies
FILTLONG=/space/sharedbin/bin/filtlong;
CUTADAPT=/space/users/smk/bin/cutadapt;
PORECHOP=/space/users/smk/bin/porechop; 
RACON=/space/sharedbin/bin/racon;
MINIASM=/space/sharedbin/bin/miniasm;
MINIMAP2=/space/users/smk/Software/minimap2-2.5/minimap2-2.5_x64-linux/minimap2;
SAMTOOLS=/space/users/smk/Software/SAMtools/bin/samtools;
PARALLEL=/usr/bin/parallel;
PRODIGAL=/space/sharedbin/bin/prodigal;
HMMSEARCH=/usr/bin/hmmsearch;
ESSENTIAL=/space/users/smk/Software/databases/essential/essential.hmm;
KAIJU=/space/users/smk/Software/kaiju/bin;
KAIJU_DB=/space/users/smk/Software/kaiju/database;
BARRNAP=/space/users/smk/Software/barrnap-0.8/bin/barrnap;
MOTHUR=/space/sharedbin/bin/mothur;
MOTHUR_DB=/space/users/smk/Software/databases/mothur/silva.seed_v132.align;
MOTHUR_TAX=/space/users/smk/Software/databases/mothur/silva.seed_v132.tax;
RSCRIPT=/usr/bin/Rscript;

### Terminal Arguments ---------------------------------------------------------


# Import user arguments
while getopts ':hzd:n:m:i:pt:l:xk' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    d) DATA_DIR=$OPTARG;;
    n) NP_ASMB=`cat $OPTARG`;;
    m) NP_COV=`cat $OPTARG`;;
    i) ILM_COV=`cat $OPTARG`;;
    f) PF=YES;;
    t) THREADS=$OPTARG;;
    l) NP_MINLENGTH=$OPTARG;;
    x) DATA_DIR=data;NP_ASMB=`cat np_asmb.txt`;NP_COV=`cat np_cov.txt`; ILM_COV=`cat ilm_cov.txt`;THREADS=40;;
    k) KEEP="YES";;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${DATA_DIR+x} ]; then echo "-d $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${NP_ASMB+x} ]; then echo "-n $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${NP_COV+x} ]; then echo "-m $MISSING"; echo "$USAGE"; exit 1; fi; 
if [ -z ${ILM_COV+x} ]; then echo "-i $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${PF+x} ]; then PF=NO; fi;
if [ -z ${THREADS+x} ]; then THREADS=10; fi;
if [ -z ${NP_MINLENGTH+x} ]; then NP_MINLENGTH=8000; fi;
if [ -z ${KEEP+x} ]; then KEEP="NO"; fi;

### Data names
NP_ALL=`printf "$NP_COV""\n""$NP_ASMB""\n" | sort | uniq`
NP_COV_NAME=`echo "$NP_COV" | sed 's/\..*$//g'`
NP_ASMB_NAME=`echo "$NP_ASMB" | sed 's/\..*$//g'`
ILM_COV_NAME=`echo "$ILM_COV" | sed 's/[12]\..*$//g' | sort | uniq`

################################################################################
### Workflow -------------------------------------------------------------------
################################################################################

### Read trimming and filtering-------------------------------------------------
if [ ! -d "trimmed_data" ]; then
mkdir -p trimmed_data

# Untrimmed data processing
if [ "$PF" = "NO" ]; then

# Nanopore data
for DATA_FILE in $NP_ALL
do
  DATA_NAME=`echo "$DATA_FILE" | sed 's/\..*$//g'`
  PORETHREADS=$(($THREADS*100000000/$(stat --printf="%s" $DATA_DIR/$DATA_FILE)))
  if [ "$PORETHREADS" -eq "0" ]; then PORETHREADS=1; fi;
  cat $DATA_DIR/$DATA_FILE | $PARALLEL --progress -j $THREADS \
  -L 4 --round-robin --pipe "cat > trimmed_data/{#}.tmp; \
  $FILTLONG --min_length $NP_MINLENGTH --min_mean_q 85 trimmed_data/{#}.tmp \
  > trimmed_data/{#}_filt.tmp; \
  $PORECHOP -i trimmed_data/{#}_filt.tmp -o trimmed_data/{#}_trim.tmp \
  --threads $PORETHREADS --min_split_read_size $NP_MINLENGTH --check_reads 1000;"
  cat trimmed_data/*_trim.tmp > trimmed_data/${DATA_NAME}_trim.fq
  rm trimmed_data/*.tmp
done

# Illumina data
for DATA_NAME in $ILM_COV_NAME
do
  $CUTADAPT -a $NEX_ADP1 -a $TRU_ADP1 -A $NEX_ADP2 \
  -A $TRU_ADP2 -j $THREADS -m 100 -q 20 \
  -o trimmed_data/${DATA_NAME}1_trim.fq \
  --paired-output trimmed_data/${DATA_NAME}2_trim.fq \
  $DATA_DIR/${DATA_NAME}1.* $DATA_DIR/${DATA_NAME}2.*
done
fi

# Trimmed data processing
if [ "$PF" = "YES" ]; then
  for DATA_FILE in $DATA_DIR/*
  do
    if [ ${DATA_FILE: -8} != "_filt.fq" ]; then
      RENAME=`basename $DATA_FILE _trim.fq`
      cp $DATA_FILE trimmed_data/${RENAME}_trim.fq
    fi    
  done
fi

# Concatenate assembly data
for DATA_NAME in $NP_ASMB_NAME
do
  cat trimmed_data/${DATA_NAME}_trim.fq >> trimmed_data/np_all_trim.fq
done
fi

### Metagenome assembly --------------------------------------------------------
if [ ! -d "metagenome_assembly" ]; then
mkdir -p metagenome_assembly

# Assembly
$MINIMAP2 -t $THREADS -x ava-ont trimmed_data/np_all_trim.fq \
trimmed_data/np_all_trim.fq > metagenome_assembly/read_minimap.paf
$MINIASM -f trimmed_data/np_all_trim.fq metagenome_assembly/read_minimap.paf \
> metagenome_assembly/assembly_miniasm.gfa
awk '/^S/{print ">"$2"\n"$3}' metagenome_assembly/assembly_miniasm.gfa |\
awk '/^>/{print ">" ++i; next}{print}' > metagenome_assembly/assembly_miniasm.fa

# Polishing
$MINIMAP2 -x map-ont -t $THREADS metagenome_assembly/assembly_miniasm.fa \
trimmed_data/np_all_trim.fq > metagenome_assembly/assembly_minimap.paf
racon -t $THREADS trimmed_data/np_all_trim.fq \
metagenome_assembly/assembly_minimap.paf metagenome_assembly/assembly_miniasm.fa \
metagenome_assembly/assembly_racon.fa
sed -i 's/^>Consensus_/>/' metagenome_assembly/assembly_racon.fa

# Circular info
echo "scaffold,circular" > metagenome_assembly/circular.csv
awk '/^S/{ print  substr( $2, 4, 6)","substr( $2, length($2), 1)}' \
metagenome_assembly/assembly_miniasm.gfa | sed -e 's/^0*//' >> metagenome_assembly/circular.csv
fi

### Metagenome mapping ---------------------------------------------------
if [ ! -d "metagenome_mapping" ]; then
mkdir -p metagenome_mapping

# Nanopore data
for DATA_FILE in $NP_COV_NAME
do
  $MINIMAP2 -ax map-ont -t $THREADS metagenome_assembly/assembly_racon.fa \
  trimmed_data/${DATA_FILE}_trim.fq | \
  $SAMTOOLS view -@ $THREADS -q 10 -Sb -F 0x104 - |\
  $SAMTOOLS sort -@ $THREADS - > metagenome_mapping/${DATA_FILE}_cov.bam
  # $SAMTOOLS view: -q 5 = min map quality; -F 0x104 = exclude unmapped and 
  # secondary hits
done

# Illumina
for DATA_FILE in $ILM_COV_NAME
do
  $MINIMAP2 -ax sr -t $THREADS metagenome_assembly/assembly_racon.fa \
  trimmed_data/${DATA_FILE}1_trim.fq trimmed_data/${DATA_FILE}2_trim.fq | \
  $SAMTOOLS view -@ $THREADS -q 10 -Sb -F 0x104 - |\
  $SAMTOOLS sort -@ $THREADS - > metagenome_mapping/${DATA_FILE}_cov.bam
done

# Read coverage and stdev
for BAM in metagenome_mapping/*_cov.bam
do
  BAM_NAME=`basename $BAM _cov.bam`
  echo "scaffold,coverage,stdev" > metagenome_mapping/${BAM_NAME}_cov.csv
  $SAMTOOLS depth -aa $BAM |\
  awk -F  "\t" '{a[$1] += $3; b[$1]++; c[$1]+=$3*$3}\
  END{OFS = ","; for (i in a) print i, a[i]/b[i], sqrt(c[i]/b[i] - (a[i]/b[i])^2)}' - |\
  sort -t, -n -k1,1 >> metagenome_mapping/${BAM_NAME}_cov.csv
  # $SAMTOOLS depth: -aa = keep contigs with 0 coverage
done

# Assembly connections
echo "scaffold1,scaffold2,connections" > metagenome_assembly/asmb_link.csv
cat metagenome_assembly/assembly_miniasm.gfa | \
grep -P "^L\tutg[0-9]{1,6}l" | \
cut -f2,4 | sed -e 's/utg0*//g' -e 's/l\t/,/g' -e 's/l$//g' -e 's/$/,5/' >> metagenome_assembly/asmb_link.csv


# Nanopore connections NB: CURRENTLY NOT WORKING
#for BAM in $NP_COV_NAME
#do
#  echo "scaffold1,scaffold2,connections" > metagenome_mapping/${BAM}_link.csv
#  CON=`$SAMTOOLS view --threads $THREADS -q 10 -F 0x104 metagenome_mapping/${BAM}_cov.bam |\
#  grep -v "^@" | cut -f1,3 | sort -u -k1,2`
#  join -j 1 -o 1.1,1.2,2.2  <(echo "$CON") <(echo "$CON") |\
#  awk '($2 < $3) {print $1,$2,$3}; ($2 > $3) {print $1,$3,$2}' |\
#  sort -u -k1,3 | cut -d" " -f2,3 | sort | uniq -c | awk '{OFS=","; print $2,$3,$1}' \
#  >> metagenome_mapping/${BAM}_link.csv
#done

# Illumina connections
for BAM in $ILM_COV_NAME
do
  echo "scaffold1,scaffold2,connections" > metagenome_mapping/${BAM}_link.csv
  CON=`$SAMTOOLS view --threads $THREADS -q 10 -F 0x104 metagenome_mapping/${BAM}_cov.bam |\
  grep -v "^@" | cut -f1,3 | sort -u -k1,2`
  join -j 1 -o 1.1,1.2,2.2  <(echo "$CON") <(echo "$CON") |\
  awk '($2 < $3) {print $1,$2,$3}; ($2 > $3) {print $1,$3,$2}' |\
  sort -u -k1,3 | cut -d" " -f2,3 | sort | uniq -c | awk '{OFS=","; print $2,$3,$1}' \
  >> metagenome_mapping/${BAM}_link.csv
done
fi

### Metagenome annotation -----------------------------------------------------
if [ ! -d "metagenome_annotation" ]; then
mkdir -p metagenome_annotation
SCF=`grep "^>" metagenome_assembly/assembly_racon.fa | tr -d ">" | sort`

# Prodigal gene prediction
cat metagenome_assembly/assembly_racon.fa | $PARALLEL --progress -j $THREADS \
--recstart '>' --pipe "cat > metagenome_annotation/{#}.tmp; $PRODIGAL \
-a metagenome_annotation/{#}_orfs.tmp -i metagenome_annotation/{#}.tmp \
-m -o /dev/null -p meta -q"
cat metagenome_annotation/*_orfs.tmp > metagenome_annotation/orfs.faa
rm metagenome_annotation/*.tmp

# Essential genes
$HMMSEARCH --tblout metagenome_annotation/hmm.orfs.txt --cut_tc --notextw \
$ESSENTIAL metagenome_annotation/orfs.faa > /dev/null
echo "scaffold,orf,hmm.id" > metagenome_annotation/essential.csv
tail -n+4  metagenome_annotation/hmm.orfs.txt | sed 's/ * / /g' | \
cut -f1,4 -d " " | sed 's/_/ /' | tr " " "," >> metagenome_annotation/essential.csv

# Kaiju protein taxonomic classification
$KAIJU/kaiju -p -z $THREADS -t $KAIJU_DB/nodes.dmp -f $KAIJU_DB/kaiju_db.fmi \
-i metagenome_annotation/orfs.faa -o metagenome_annotation/kaiju.out
$KAIJU/addTaxonNames -u -r phylum -t $KAIJU_DB/nodes.dmp -n $KAIJU_DB/names.dmp \
-i metagenome_annotation/kaiju.out -o metagenome_annotation/kaiju.names.out

# Majority vote contig classification
echo "scaffold,phylum" > metagenome_annotation/tax.csv
cat metagenome_annotation/kaiju.names.out | \
sed -e 's/_/\t/' -e '/NA;/d' -e 's/; //'  | \
cut -f2,5 | \
awk -F "\t" '{a[$1","$2]++} END{OFS = ","; for (i in a) print i, a[i]}' - | \
awk -F "," '{if (a[$1]<$3){a[$1]=$1; b[$1]=$2; c[$1]=$3}; d[$1]+=$3} \
END{OFS = ","; for (i in a){if (c[i] >= 2 && c[i]/d[i] > 0.6) print a[i], b[i], c[i] }}' - |\
sort -t, -k1,1 | \
join -a 1 -e none -t, -j 1 -o 1.1,2.2 <(echo "$SCF") - | sort -n -t, -k1,1 \
>> metagenome_annotation/tax.csv

# 16S rRNA gene extraction
$BARRNAP metagenome_assembly/assembly_racon.fa --reject 0.3 --threads $THREADS \
--kingdom bac --quiet > metagenome_annotation/rRNA_search.txt
grep "16S_rRNA" metagenome_annotation/rRNA_search.txt | cut -f1,4,5 | \
sed -e '/^#/d' -e 's/Name=//' | sort -u > metagenome_annotation/ssu_gene.txt

$SAMTOOLS faidx metagenome_assembly/assembly_racon.fa
for GENE in "$(cat metagenome_annotation/ssu_gene.txt)"
do
 REGION=`echo "$GENE" | sed -e 's/\t/:/' -e 's/\t/-/'`
 $SAMTOOLS faidx metagenome_assembly/assembly_racon.fa $REGION \
 >> metagenome_annotation/ssu_gene.fa
done

$MOTHUR "#classify.seqs(fasta=metagenome_annotation/ssu_gene.fa,
reference=$MOTHUR_DB, taxonomy=$MOTHUR_TAX,
processors=$THREADS, outputdir=./metagenome_annotation)"

echo "scaffold,ssu" > metagenome_annotation/ssu_tax.csv
cat metagenome_annotation/ssu_gene.seed_v132.wang.taxonomy | \
sed 's/_[0-9-]*\t/,/' >> metagenome_annotation/ssu_tax.csv

echo "scaffold,ssu_count" > metagenome_annotation/ssu_count.csv
cat metagenome_annotation/ssu_tax.csv | tail -n +2 | cut -d"," -f1 | \
sort | uniq -c | sed -e "s/^ *//" | awk '{ print $2 "," $1}' | \
join -a 1 -e 0 -t, -j 1 -o 1.1,2.2 <(echo "$SCF") - | sort -t, -k1,1 -n \
>> metagenome_annotation/ssu_count.csv
fi


### Prepare for binning --------------------------------------------------------
if [ ! -d "binning" ]; then

mkdir -p binning
mkdir -p binning/data
mkdir -p binning/bins

cp metagenome_assembly/assembly_racon.fa binning/data/
cp metagenome_assembly/circular.csv binning/data/
cp metagenome_assembly/asmb_link.csv binning/data/
cp metagenome_mapping/*_cov.csv binning/data/
cp metagenome_mapping/*_link.csv binning/data/
cp metagenome_annotation/tax.csv binning/data/
cp metagenome_annotation/essential.csv binning/data/
cp metagenome_annotation/ssu_tax.csv binning/data/
cp metagenome_annotation/ssu_count.csv binning/data/


echo "#mmlong data loading:
setwd(\"./binning\")

require(\"mmgenome\", quietly = T)
require(\"tidyverse\", quietly = T)
options(scipen = 8)

assembly <- readDNAStringSet(\"data/assembly_racon.fa\", format = \"fasta\")

cov_list <- list.files(path = \"data/\", pattern = \".*_cov.csv\", full.names = T)
for (i in cov_list) {
  assign(gsub(\".csv\", \"\", basename(i)), read_csv(i, T)[,1:2])
}

link_list <- list.files(path = \"data/\", pattern = \".*_link.csv\", full.names = T)
for (i in link_list) {
  assign(gsub(\".csv\", \"\", basename(i)), read_csv(i, T))
}

circ <- read_csv(\"data/circular.csv\", T)

tax <- read_csv(\"data/tax.csv\", T)

ess <- read_csv(\"data/essential.csv\", T, comment = \"#\")

ssu <- read_csv(\"data/ssu_count.csv\", T)

d <- mmload(assembly = assembly, pca = T, BH_SNE = F, threads = 8,
            coverage = gsub(\".csv\", \"\", basename(cov_list)),
            ess = ess)
d\$scaffolds\$tax <- tax\$phylum
d\$scaffolds\$ssu <- ssu\$ssu_count
d\$scaffolds\$circ <- circ\$circular

rm(list = grep(c(\"cov|tax|link_list|first_time|ess|circ|ssu\"), ls(), value = T))
rm(i)
save.image(file=\"data.RData\")
" > binning/load.R

$RSCRIPT binning/load.R

echo "# mmlong binning template

\`\`\`{r}
## Import packages and data
library(\"mmgenome\")
library(grid)
options(scipen = 8)
mmimport(\"Load_data.Rmd\")
\`\`\`

\`\`\`{r}
## Overview plot
mmplot(data = d, x = \"np_cov\", y = \"ilm_cov\", color = \"gc\") +
  theme(legend.position = \"right\")
\`\`\`
" > binning/analysis.Rmd
fi

### Cleanup --------------------------------------------------------------------
if [ "$KEEP" = "NO" ]; then
  rm -rf ./metagenome_annotation
  rm -rf ./metagenome_assembly
  rm -f ./metagenome_mapping/*.csv
fi


### Logfile --------------------------------------------------------------------
LOG_NAME="mmlong_metagenome_log_$(date +%s).txt"

echo "Input arguments:
mmlong metagenome workflow script version: $VERSION
Sequencing data folder: $DATA_DIR
Nanopore data list for assembly: $NP_ASMB
Nanopore data list for coverage estimates: $NP_COV
Illumina data list for coverage estimates: $ILM_COV
Number of threads: $THREADS
Nanopore min read length: $NP_MINLENGTH" >> $LOG_NAME

################################################################################
### Exit -----------------------------------------------------------------------
################################################################################
exit 0

################################################################################
### Testing --------------------------------------------------------------------
################################################################################
NP_ASMB=`cat np_asmb.txt`;
NP_COV=`cat np_cov.txt`;
ILM_COV=`cat ilm_cov.txt`;
ILM_ADP1=CTGTCTCTTATACACATCTGACGCTGCCGACGA
ILM_ADP2=CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
FILTLONG=/space/sharedbin/bin/filtlong;
CUTADAPT=/space/users/smk/bin/cutadapt;
PORECHOP=/space/users/smk/bin/porechop; 
RACON=/space/sharedbin/bin/racon;
MINIASM=/space/sharedbin/bin/miniasm;
MINIMAP2=/space/users/smk/Software/minimap2-2.5/minimap2-2.5_x64-linux/minimap2;
SAMTOOLS=/space/users/smk/Software/SAMtools/bin/samtools;
PARALLEL=/usr/bin/parallel;
PRODIGAL=/space/sharedbin/bin/prodigal;
HMMSEARCH=/usr/bin/hmmsearch;
ESSENTIAL=/space/users/smk/Software/databases/essential/essential.hmm;
KAIJU=/space/users/smk/Software/kaiju/bin;
KAIJU_DB=/space/users/smk/Software/kaiju/database;
BARRNAP=/space/users/smk/Software/barrnap-0.8/bin/barrnap;
MOTHUR=/space/sharedbin/bin/mothur;
MOTHUR_DB=/space/users/smk/Software/databases/mothur/silva.seed_v132.align;
MOTHUR_TAX=/space/users/smk/Software/databases/mothur/silva.seed_v132.tax;
RSTUDIO=/usr/bin/rstudio;
NP_ALL=`printf "$NP_COV""\n""$NP_ASMB""\n" | sort |uniq`;
NP_COV_NAME=`echo "$NP_COV" | sed 's/\..*$//g'`;
NP_ASMB_NAME=`echo "$NP_ASMB" | sed 's/\..*$//g'`;
ILM_COV_NAME=`echo "$ILM_COV" | sed 's/[12]\..*$//g' | sort | uniq`;
THREADS=40;
NP_MINLENGTH=8000;
DATA_DIR=../data;


