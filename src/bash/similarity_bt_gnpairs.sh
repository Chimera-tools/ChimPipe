#!/bin/bash

# Difference with previous version is the fact that exonic length of each transcript is not provided
# and that final file is not gzipped (in order to be used by other programs)

# usage
#######
# similarity_bt_gnpairs.sh annot genome_GEM threads

# example
#########
# time similarity_bt_gnpairs.sh gen10.long.exon.gtf hg19.gem 4

# Notes
#######
# - Made for using on a 64 bit linux architecture


# Programs
##########
EXTRACT_TRANSC_SEQS=/users/rg/sdjebali/bin/gtf2fasta.sh
MAKEDB=/users/rg/sdjebali/bin/ncbi-blast-2.2.29+/bin/makeblastdb
BLAST=/users/rg/sdjebali/bin/ncbi-blast-2.2.29+/bin/blastn


# In case the user does not provide the two necesary input files
################################################################
if [ ! -n "$1" ] || [ ! -n "$2" ]
then
    echo "" >&2
    echo Usage: similarity_bt_gnpairs.sh annot genome_GEM threads >&2
    echo "" >&2
    echo Example: similarity_bt_gnpairs.sh gen10.long.exon.gtf hg19.gem 4 >&2
    echo "" >&2
    exit 1
else
	annot=$1
	index=$2
fi

if [ ! -n "$3" ]
then
    threads=1;
else
    threads=$3;
fi

# Variable from input
#####################
b=`basename $annot`
b2tmp=${b%.gtf}
b2=${b2tmp%.gff}


start=$(date +%s)

# 1. extract the cdna sequence of each transcript in the annotation
###################################################################
echo I am extracting the cdna sequence of each transcript \in the annotation >&2
$EXTRACT_TRANSC_SEQS $annot $index 
echo done >&2

# 2. make BLAST DB out of the transcript sequences
##################################################
echo I am making a BLAST database out of the transcript sequences >&2
$MAKEDB -in $b2\_tr.fasta -input_type 'fasta' -dbtype 'nucl' -parse_seqids 
echo done >&2

# 3. run Blast on all against all to detect local similarity between transcripts (long)
#######################################################################################
echo I am running Blast on all against all to detect local similarity between transcripts >&2
$BLAST -query $b2\_tr.fasta -db $b2\_tr.fasta -lcase_masking -dust 'no' -ungapped -outfmt '7' -num_threads $threads | gzip > $b2\_tr_vs_$b2\_tr.blastout.tsv.gz
echo done >&2

# 4. get a gene pair file with % similarity, alignment length and other information (10 minutes)
##################################################################################################
echo I am making a gene pair file with % similarity, alignment length and other information >&2
zcat $b2\_tr_vs_$b2\_tr.blastout.tsv.gz | awk -v fileRef=$annot 'BEGIN{while (getline < fileRef >0){if($3=="exon"){split($10,a,"\""); split($12,b,"\""); gene[b[2]]=a[2]}}} $1!~/#/{gp=((gene[$1]<=gene[$2]) ? (gene[$1]"-"gene[$2]) : (gene[$2]"-"gene[$1])); if((sim[gp]=="")||(sim[gp]<=$3)){if((lgal[gp]=="")||(lgal[gp]<$4)){sim[gp]=$3; lgal[gp]=$4; tr[gp]=((gene[$1]<=gene[$2]) ? ($1"-"$2) : ($2"-"$1))}}} END{for(gp in sim){split(gp,a,"-"); split(tr[gp],b,"-"); print a[1], a[2], sim[gp], lgal[gp], b[1], b[2]}}' | awk '$1!=$2{print $1, $2, $3, $4, $5, $6}' | sort | uniq > $b2"_"gene1_gene2_alphaorder_pcentsim_lgalign_trpair.txt
echo done >&2

end=$(date +%s)
echo "Completed in $(echo "($start-$end)/60" | bc -l | xargs printf "%.2f\n") min" >&2

