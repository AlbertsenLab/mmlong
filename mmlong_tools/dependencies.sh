#!/bin/bash
# mmlong dependencies
# By Søren Karst
# Version
VERSION=0.1.1

### Paths to dependencies

# General
FILTLONG=/space/sharedbin/bin/filtlong;
RACON=/space/sharedbin/bin/racon;
MINIASM=/space/sharedbin/bin/miniasm;
MINIMAP2=/space/users/smk/Software/minimap2-2.5/minimap2-2.5_x64-linux/minimap2;
SAMTOOLS=/space/users/smk/Software/SAMtools/bin/samtools;
PARALLEL=/usr/bin/parallel;

# Metaflow
CUTADAPT=/space/users/smk/bin/cutadapt;
PORECHOP=/space/users/smk/bin/porechop; 
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

# Reassembly
SPADES=/space/users/smk/Software/Spades3.11/SPAdes-3.11.1-Linux/bin/spades.py;
UNICYCLER=/space/users/smk/bin/unicycler;

# Stats
QUAST=/space/sharedbin/bin/quast;

# Legacy
#MAUVE=/space/users/smk/Software/mauve_snapshot_2015-02-13/Mauve.jar;
#PMAUVE=/space/users/smk/Software/mauve_snapshot_2015-02-13/linux-x64/progressiveMauve;
