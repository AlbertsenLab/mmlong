#!/bin/bash
# mmlong dependencies
# By Søren Karst
# Version
VERSION=0.1.4

### Paths to dependencies
# mmlong
MMLONG_LINKS=/space/users/smk/Desktop/2017_long_read/pipeline/mmlong/mmlong_tools/mmlong_links;
MMLONG_READCOVERAGE=/space/users/smk/Desktop/2017_long_read/pipeline/mmlong/mmlong_tools/mmlong_readcoverage;
MMLONG_MINIASM_RACON=/space/users/smk/Desktop/2017_long_read/pipeline/mmlong/mmlong_tools/mmlong_miniasm-racon;

# General
FILTLONG=/space/sharedbin/bin/filtlong;
RACON=/space/users/smk/Software/racon/build/bin/racon;
MINIASM=/space/sharedbin/bin/miniasm;
MINIMAP2=//space/users/smk/Software/minimap2-2.11/minimap2;
SAMTOOLS=/space/users/smk/Software/SAMtools/bin/samtools;
PARALLEL=/usr/bin/parallel;
SPADES=/space/users/smk/Software/Spades3.11/SPAdes-3.11.1-Linux/bin/spades.py;

# Metaflow
CUTADAPT=/usr/local/bin/cutadapt;
PORECHOP=/space/users/smk/bin/porechop;
PRODIGAL=/space/sharedbin/bin/prodigal;
FRAGGENESCAN=/space/users/smk/Software/FragGeneScan1.30/run_FragGeneScan.pl;
HMMSEARCH=/usr/bin/hmmsearch;
ESSENTIAL=/space/users/smk/Desktop/2017_long_read/pipeline/mmlong/mmlong_tools/databases/essential.hmm;
KAIJU=/space/users/smk/Software/kaiju/bin;
KAIJU_DB=/space/users/smk/Software/kaiju/database;
BARRNAP=/space/users/smk/Software/barrnap-0.8/bin/barrnap;
MOTHUR=/space/sharedbin/bin/mothur;
MOTHUR_DB=/space/users/smk/Desktop/2017_long_read/pipeline/mmlong/mmlong_tools/databases/silva.seed_v132.align;
MOTHUR_TAX=/space/users/smk/Desktop/2017_long_read/pipeline/mmlong/mmlong_tools/databases/silva.seed_v132.tax;
RSCRIPT=/usr/bin/Rscript;
METABAT2_COV=/space/sharedbin/bin/jgi_summarize_bam_contig_depths;
METABAT2=/space/sharedbin/bin/metabat2;
MEGAHIT=/space/users/smk/bin/megahit;

# Reassembly
UNICYCLER=/space/users/smk/bin/unicycler;

# Stats
QUAST=/space/sharedbin/bin/quast;

# Legacy
#MAUVE=/space/users/smk/Software/mauve_snapshot_2015-02-13/Mauve.jar;
#PMAUVE=/space/users/smk/Software/mauve_snapshot_2015-02-13/linux-x64/progressiveMauve;
