#!/bin/bash

<<authors
*****************************************************************************
	
	ChimPE.sh
	
	This file is part of the ChimPipe pipeline 

	Copyright (c) 2014 Bernardo Rodríguez-Martín 
					   Emilio Palumbo 
					   Sarah Djebali 
	
	Computational Biology of RNA Processing group
	Department of Bioinformatics and Genomics
	Centre for Genomic Regulation (CRG)
					   
	Github repository - https://github.com/Chimera-tools/ChimPipe
	
	Documentation - https://chimpipe.readthedocs.org/

	Contact - chimpipe.pipeline@gmail.com
	
	Licenced under the GNU General Public License 3.0 license.
******************************************************************************
authors

# Description
##############

# usage
#######
# bash ChimPE.sh alignments.bam genome_index.gem annot.gtf normalJunctions_ChimSplice readDirectionality outDir

# Input
########
# 1) BAM file. Mandatory
# 2) GEM indexed genome. Mandatory
# 3) reference annotation. Mandatory
# 4) normal splice junctions produced by ChimSplice. Mandatory 
# 5) readDirectionality (MATE1_SENSE|MATE2_SENSE|UNSTRANDED). Default: UNSTRANDED
# 6) Output directory. Default: current working directory

# Output
###########
# 1) discordant_readPairs.txt
# SRR201779.4705817_PATHBIO-SOLEXA2_30TUEAAXX:3:73:27:1255_length=53#0 100 ENSG00000157554.14:chr21_39817327_39817544,ENSG00000157554.14:chr21_39817327_39817544,ENSG00000157554.14:chr21_39817327_39817544 ENSG00000157554.14 ERG protein_coding 98.1132 ENSG00000184012.7:chr21_42880008_42880086,ENSG00000184012.7:chr21_42880008_42880085 ENSG00000184012.7 TMPRSS2 protein_coding


# will exit if there is an error or in a pipe
set -e -o pipefail

# In case the user does not provide any input file
###################################################
if [[ ! -e "$1" ]] || [[ ! -e "$2" ]] || [[ ! -e "$3" ]]
then
    echo "" >&2
    echo "*** ChimPE V0.9.0 ***" >&2
    echo "" >&2
    echo "Usage:    bash ChimPE.sh alignments.bam genome_index.gem annot.gtf normalJunctions_ChimSplice readDirectionality outDir" >&2
    echo "" >&2
    echo "Example:  bash ChimPE.sh vcap.bam Homo_sapiens.GRCh37.chromosomes.chr.gem gencode.v7.annotation_exons.gtf normal_spliceJunctions.txt UNSTRANDED ~brodriguez/projects/VcaP" >&2
    echo "" >&2
    echo "" >&2
    echo "" >&2
    echo "" >&2
    echo "Input:" >&2
    echo "1) BAM file. Mandatory" >&2
    echo "2) GEM indexed genome. Mandatory" >&2
    echo "3) reference annotation. Mandatory" >&2
    echo "4) normal splice junctions produced by ChimSplice. Mandatory" >&2
    echo "5) readDirectionality (MATE1_SENSE|MATE2_SENSE|UNSTRANDED). Default: UNSTRANDED"  >&2
    echo "6) Output directory. Default: current working directory"  >&2
    echo "" >&2
    echo "Output:" >&2 
    echo "1) discordant_readPairs.txt" >&2
    echo "" >&2
    exit -1
fi

# In case the user does not provide any indexed genome, annotation 
##################################################################
# file or output dir or strandedness we provide default values
##############################################################
bam=$1;
genome=$2;
annot=$3;
normalJunc=$4

if [ ! -n "$5" ]
then
	readDirectionality="UNSTRANDED"
	outDir=.
else
	readDirectionality=$5
	if [ ! -n "$6" ]
	then
		outDir=.
	else
		outDir=$6
	fi
fi	


# Directories 
#############
## Set root directory
path="`dirname \"$0\"`"              # relative path
rootDir="`( cd \"$path\" && pwd )`"  # absolute path

if [ -z "$rootDir" ] ; 
then
  # error; for some reason, the path is not accessible
  # to the script
  log "Path not accessible to the script\n" "ERROR" 
  exit 1  # fail
fi

## Set bin and awk scripts directories
binDir=$rootDir/../../bin
awkDir=$rootDir/../awk

# Programs and scripts
########################

## Bin
GEMINFO=$binDir/gemtools-1.7.1-i3/gem-info

## Awk
BED2GFF=$awkDir/bed2gff.awk
GFF_CORRECT_STRAND=$awkDir/gffCorrectStrand.awk
GFF2GFF=$awkDir/gff2gff.awk
PERC_OVERLAP=$awkDir/compute_mappedReadPercOverlap2Exon.awk
SELECT_OVERLAPPING_EXONS=$awkDir/select_overlappingExons_contiguousMapping.awk
CLASSIFY_PE=$awkDir/classifyPairedEnds.awk
DETECT_DISCORDANT=$awkDir/detect_discordantPEs_1mateSplitmap_1mateMapContiguously.awk


# 0) GENERATE GENOME FILE # 
###########################
# - $outDir/chromosomes_length.txt

$GEMINFO $genome | awk -v OFS='\t' '$1 ~ /^#/ && ($2 !~ /M/) && ($2 !~ /Mt/) && ($2 !~ /MT/) && ($3 ~ /F/) {split($5,lengths,"@"); chrLengths[$2]=lengths[2];}END{for (chr in chrLengths){print chr, chrLengths[chr]+1;}}' | tr -d \'\" | sort -V -k1,1 > $outDir/chromosomes_length.txt


#####################################################
# 1) MAKE GFF CONTAINING THE ALIGNMENT COORDINATES  #
#    FOR THE CONTIGUOUSLY MAPPED READS              #
#####################################################
# - $outDir/contiguousAlignments.gff.gz
# NOTE: for now it only considers unique mappings. Multimapped reads are discarded. 

echo "1. Make GFF containg the alingment coordinates for the contiguously mapped reads" >&2

## Only unique mappings considered ($5==1)
samtools view -bF 4 $bam | bedtools bamtobed -bed12 -tag NH -i stdin | awk '($5==1) && ($10==1) && ($1 !~ /M/) && ($1 !~ /Mt/) && ($1 !~ /MT/)' | awk -f $BED2GFF | awk -v OFS='\t' -v readDirectionality=$readDirectionality -f $GFF_CORRECT_STRAND | awk -f $GFF2GFF | gzip > $outDir/contiguousAlignments.gff.gz
	

#########################################################
# 2) INTERSECT READ ALIGNMENTS WITH THE ANNOTATED EXONS #
#########################################################
# - $outDir/contiguousAlignments_intersected.txt.gz

echo "2. Intersect read alignments with the annotated exons" >&2

bedtools intersect -wao -a $outDir/contiguousAlignments.gff.gz -b $annot -sorted -g $outDir/chromosomes_length.txt | gzip > $outDir/contiguousAlignments_intersected.txt.gz

###########################################################
# 3) FOR EACH INTERSECTION ADD THE PERCENTAGE OF OVERLAP  #
# BETWEEN THE MAPPED READ AND THE OVERLAPPING EXON        #
###########################################################
# - $outDir/spliceJunctions_2parts_beg_end_intersected_percOverlap_distExonSS.txt.gz

echo "3. For each intersection add the percentage of overlap between the mapped read and the overlapping exon " >&2

awk -f $PERC_OVERLAP <( gzip -dc $outDir/contiguousAlignments_intersected.txt.gz ) | gzip > $outDir/contiguousAlignments_intersected_percOverlap.txt.gz


#########################################################################
# 4) FOR EACH MAPPED READ SELECT THE OVERLAPPING EXONS WHICH BETTER FIT #
#########################################################################
# - $outDir/contiguousAlignments_intersected_percOverlap_exonsBetterFit.txt.gz

echo "4. For each mapped read select those overlaping exons which maximize the percentage of overlap" >&2

awk -f $SELECT_OVERLAPPING_EXONS <( gzip -dc $outDir/contiguousAlignments_intersected_percOverlap.txt.gz ) | gzip > $outDir/contiguousAlignments_intersected_percOverlap_exonsBetterFit.txt.gz


#################################################################################################
# 5) CLASSIFY READ PAIRS INTO: 																	#
#     A) CONCORDANT:  BOTH MATES MAPPING IN ONLY ONE GENE										#
#     B) DISCORDANT:  BOTH MATES MAPPING IN DIFFERENT GENES										#
#     C) UNANNOTATED: BOTH MATES MAP AND AT LEAST ONE OF THEM DO NOT OVERLAP ANY ANNOTATED GENE #
#     D) UNPAIRED:    ONLY ONE PAIR MAPPING														#
#################################################################################################
# - $outDir/discordantPE.txt

echo "5. Classify read pairs in concordant, discordant, unannotated or unpaired" >&2

awk -v OFS='\t' -f $CLASSIFY_PE <( gzip -dc $outDir/contiguousAlignments_intersected_percOverlap_exonsBetterFit.txt.gz ) > $outDir/readPairs_classified.txt

### Make a different file for each type of read pair:

## A) CONCORDANT
awk '($2=="CONCORDANT") || (NR==1)' $outDir/readPairs_classified.txt > $outDir/concordant_contiguousMapped_readPairs.txt

## B) DISCORDANT
awk '($2=="DISCORDANT") || (NR==1)' $outDir/readPairs_classified.txt > $outDir/discordant_contiguousMapped_readPairs.txt

## C) UNANNOTATED
awk '($2=="UNANNOTATED") || (NR==1)' $outDir/readPairs_classified.txt > $outDir/unannotated_contiguousMapped_readPairs.txt

## D) UNPAIRED
awk '($2=="UNPAIRED") || (NR==1)' $outDir/readPairs_classified.txt > $outDir/unpaired_contiguousMapped_reads.txt


#####################################################################
# 6) FIND READ PAIRS WHERE ONE OF THE MATES MAPS CONTIGUOUSLY AND    #
# THE OTHER ONE SPANS A SPLICE JUNCTION (SPLIT-MAPS). CLASSIFY THEM #
#####################################################################

echo "6. Detect additional discordant pairs where one maps contiguously and the other one split-map" >&2

awk -v OFS='\t' -v fileRef=$normalJunc -f $DETECT_DISCORDANT $outDir/unpaired_contiguousMapped_reads.txt > $outDir/discordant_contiguousAndSplitMapped_readPairs.txt

#######################################################################
## 7) PRODUCE THE FINAL DISCORDANT OUTPUT FILE COMBINING 5.B) AND 6)  #
#######################################################################

echo "7. Produce a final discordant pairs output file" >&2

cat $outDir/discordant_contiguousMapped_readPairs.txt $outDir/discordant_contiguousAndSplitMapped_readPairs.txt > $outDir/discordant_readPairs.txt


######################
# 8) CLEANUP AND END #
######################
echo "8. Cleanup and end" >&2
rm $outDir/chromosomes_length.txt $outDir/contiguousAlignments.gff.gz $outDir/contiguousAlignments_intersected.txt.gz $outDir/contiguousAlignments_intersected_percOverlap.txt.gz $outDir/contiguousAlignments_intersected_percOverlap_exonsBetterFit.txt.gz $outDir/readPairs_classified.txt $outDir/discordant_contiguousMapped_readPairs.txt $outDir/discordant_contiguousAndSplitMapped_readPairs.txt




