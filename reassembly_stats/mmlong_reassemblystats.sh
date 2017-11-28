#!/bin/bash

#####################
#### Description ####
#####################

USAGE="$(basename "$0") [-h] [-cpath string -threads int -fileending string] -- reassembly_stats: Calculation of basic assembly stats from a folder with bins.

Requirements:
checkm (v. 1.0.7)
barrnap (v. 0.7)
gff.prokka.to.table.pl ()
parseCDS.R ()
"
###################
#### Arguments ####
###################

# PATHS
CPATH=bins/;
# Settings
THREADS=12;
FILEENDING=fasta;


##################
#### Analysis ####
##################

mkdir -p results
mkdir -p temp/rrna
mkdir -p temp/checkm_results

echo "Running CheckM"
checkm lineage_wf -x $FILEENDING -t $THREADS --tab_table $CPATH temp/checkm_results -f temp/checkm.tsv
# Add a header with names
sed -i '1d' temp/checkm.tsv
sed -i '1s/^/bin\tmarker_lineage\tn_genomes\tn_markers\tn_marker_sets\tS0\tS1\tS2\tS3\tS4\tS5plus\tCompleteness\tContamination\tStrain_heterogeneity\n/' temp/checkm.tsv

# Revisit this if using a different version to ensure that the format has not changed
echo "Parsing checkm binstats (for some reason their tsv format does not have a header)"
# Remove variable names within the lines
sed -e "s/{'Translation table': /,/g" -e "s/ 'GC std': //g" -e "s/'# ambiguous bases': //g" -e "s/'Genome size'://g" -e "s/'Longest contig'://g" -e "s/'N50 (scaffolds)'://g" -e "s/'Mean scaffold length'://g" -e "s/'# contigs'://g"  -e "s/'# scaffolds'://g" -e "s/'# predicted genes'://g" -e "s/'Longest scaffold'://g" -e "s/'GC'://g" -e "s/'N50 (contigs)'://g" -e "s/'Coding density'://g" -e "s/'Mean contig length'://g" -e "s/}//g" -e "s/\s//g" temp/checkm_results/storage/bin_stats.analyze.tsv > temp/checkm_binstats_parsed.csv
# Add a header with names
sed -i '1s/^/bin,Translation_table,GC_std,n_ambiguous_bases,Genome_size,Longest_contig,N50_scaffolds,Mean_scaffold_length,n_contigs,n_scaffolds,# predicted_genes,Longest_scaffold,GC,N50_contigs,Coding_density,Mean_contig_length\n/' temp/checkm_binstats_parsed.csv

echo "bin,medianCDSsize" > temp/CDSmedian.csv;
echo "Calculating median CDS size from checkm"
for F in $CPATH/*.$FILEENDING; do 
  N=$(basename $F .$FILEENDING) ; 
  perl gff.prokka.to.table.pl -i temp/checkm_results/bins/$N/genes.gff -o temp/$N.CDS.txt;
  Rscript parseCDS.R temp/$N.CDS.txt $N >> temp/CDSmedian.csv;
done

echo "bin,found16S" > temp/found16S.csv
echo "Use barrnap to check whether the bin has a 16S rRNA sequence"
for F in $CPATH/*.$FILEENDING; do 
  N=$(basename $F .$FILEENDING) ; 
  barrnap --threads $THREADS --kingdom bac --quiet $F > temp/rrna/$N.txt
if grep -q 16S_rRNA temp/rrna/$N.txt; then
               echo "$N,1" >> temp/found16S.csv
               echo "Found 16S in $N"
            else
               echo "$N,0" >> temp/found16S.csv
               echo "Did not find 16S in $N"
            fi
done

# Merge output to a single file
Rscript mergereassemblystats.R results/reassemblystats.csv;
echo "Done -->"

