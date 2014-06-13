#!/bin/bash


# 10/01/2013
############

# /users/rg/brodriguez/Bin/similarity_bt_gnpairs.sh 
####################################################

# usage
#######
# similarity_bt_gnpairs.sh annot chromDir 

# example
#########
# time  similarity_bt_gnpairs.sh /users/rg/projects/encode/scaling_up/whole_genome/Gencode/version10/Long/Element/gen10.long.exon.gtf /db/seq/genomes/H.sapiens/golden_path_200902/chromFa 

# Notes
#######
# - Made for using on a 64 bit linux architecture


# Programs
###########
DIVIDE_CHR=~sdjebali/bin/divide_into_chr.sh 
EXTRACT_TRANSC_SEQS=/users/rg/jlagarde/julien_utils/extract_spliced_transc_seqs.pl
FORMATDB=/users/rg/brodriguez/bin/formatdb
BLAST=/users/rg/brodriguez/bin/blastall

# In case the user does not provide any input file
###################################################
if [ ! -n "$1" -a ! -n "$2" ]
then
    echo "" >&2
    echo Usage: similarity_bt_gnpairs.sh annot chromdir >&2
    echo "" >&2
    echo Example: similarity_bt_gnpairs.sh /users/rg/projects/encode/scaling_up/whole_genome/Gencode/version10/Long/Element/gen10.long.exon.gtf /db/seq/genomes/H.sapiens/golden_path_200902/chromFa >&2
    echo "" >&2
    exit 1
else
annot=$1
chromDir=$2
fi

# 1. extract the cdna sequence of each transcript in the annotation
###################################################################
# a. divide into chr
####################

format=${annot##*.}
annotbase=`basename ${annot%"."${format}}`

###$DIVIDE_CHR $annot $format

# b. extract chr names
###################
#cat $annot | cut -f1 | sort | uniq > chr.list

# c. extract transcript sequences (long)
#################################
#cat chr.list | while read chr
#do 
##$EXTRACT_TRANSC_SEQS $annotbase"_"$chr"."$format $chromDir > $annotbase"_"$chr\_transcripts.fa
#done

# d. cat into one file and delete intermediate files
########################################
#cat chr.list | while read chr
#do 
#cat $annotbase"_"$chr\_transcripts.fa 
#done > $annotbase"_"transcripts_hg19.fa

#cat chr.list | while read chr
#do 
#rm $annotbase"_"$chr"."$format $annotbase"_"$chr\_transcripts.fa 
#done

# e. Compute the exonic length of each transcript
#########################################
#fastalength $annotbase"_"transcripts_hg19.fa | awk '{print $2, $1}' > $annotbase.tr.exlg.txt

# 2. make BLAST DB out of the transcript sequences and run Blast on all against all to detect local similarity between transcripts (long)
################################################################################################################
##$FORMATDB -V -i $annotbase"_"transcripts_hg19.fa -p F
##$BLAST -g F -a 6 -U F -p blastn -i $annotbase"_"transcripts_hg19.fa -d $annotbase"_"transcripts_hg19.fa -F F -m 9 | gzip > $annotbase"_"transcripts_hg19_vs_$annotbase"_"transcripts_hg19.blastout.tsv.gz

# 4. Get a gene pair file with % similarity, alignment length and other information (10 minutes)
#############################################################################

zcat $annotbase"_"transcripts_hg19_vs_$annotbase"_"transcripts_hg19.blastout.tsv.gz | awk -v fileRef=$annot 'BEGIN{while (getline < fileRef >0){if($3=="exon"){split($10,a,"\""); split($12,b,"\""); gene[b[2]]=a[2]}}} $1!~/#/{gp=((gene[$1]<=gene[$2]) ? (gene[$1]"-"gene[$2]) : (gene[$2]"-"gene[$1])); if((sim[gp]=="")||(sim[gp]<=$3)){if((lgal[gp]=="")||(lgal[gp]<$4)){sim[gp]=$3; lgal[gp]=$4; tr[gp]=((gene[$1]<=gene[$2]) ? ($1"-"$2) : ($2"-"$1))}}} END{for(gp in sim){split(gp,a,"-"); split(tr[gp],b,"-"); print a[1], a[2], sim[gp], lgal[gp], b[1], b[2]}}' | awk -v fileRef=$annotbase.tr.exlg.txt 'BEGIN{while (getline < fileRef >0){exlg[$1]=$2}} $1!=$2{print $1, $2, $3, $4, $5, $6, exlg[$5], exlg[$6]}' | sort | uniq | gzip > gene1_gene2_alphaorder_pcentsim_lgalign_trpair_trexoniclength.txt.gz


