#!/bin/bash

<<authors
*****************************************************************************
	
	ChimSplice.sh
	
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
# Make chimeric and normal junctions from spliced aligned reads (in BAM and GEM files), annotated them and produce a file for each type of junction
# with the number of reads supporting them, the overlapping genes and other useful info associated

# usage
#######
# bash ChimSplice.sh alignment_files_paths.txt genome_index.gem annot.gtf readDirectionality consensusSS outDir

# Input
########
# 1) file listing the absolute paths to a set of BAM or GEM alignment files (GEM produced with the rna-mapper)(one row for each path). Mandatory
# 2) GEM indexed genome. Mandatory
# 3) reference annotation. Mandatory
# 4) readDirectionality (MATE1_SENSE|MATE2_SENSE|UNSTRANDED). Default: UNSTRANDED
# 5) consensus splice sites used for splice alignments. Default='GT+AG,GC+AG,ATATC+A.,GTATC+AT'
# 6) Output directory. Default: current working directory

# Output
###########
# 1) chimJunctions_candidates.txt
# chr1_26515889_-:chr1_1677305_-	6	13	0	13	26515903	1677255	GC	AG	1	1	24838584	100	100 119	53	ENSG00000142675.13	ENSG00000215790.2	CNKSR1	SLC35E2	protein_coding	protein_coding	SRR201779.1193309_PATHBIO-SOLEXA2_30TUEAAXX:3:18:298:864_length=53#0/1,SRR201779.1193309_PATHBIO-SOLEXA2_30TUEAAXX:3:18:298:864_length=53#0/1,SRR201779.911235_PATHBIO-SOLEXA2_30TUEAAXX:3:14:750:1379_length=53#0/2,SRR201779.5697889_PATHBIO-SOLEXA2_30TUEAAXX:3:88:1730:1110_length=53#0/2,SRR201779.5697889_PATHBIO-SOLEXA2_30TUEAAXX:3:88:1730:1110_length=53#0/2,SRR201779.4370469_PATHBIO-SOLEXA2_30TUEAAXX:3:68:681:265_length=53#0/2,SRR201779.4370469_PATHBIO-SOLEXA2_30TUEAAXX:3:68:681:265_length=53#0/2,SRR201779.3242759_PATHBIO-SOLEXA2_30TUEAAXX:3:50:1449:676_length=53#0/2,SRR201779.3242759_PATHBIO-SOLEXA2_30TUEAAXX:3:50:1449:676_length=53#0/2,SRR201779.1919192_PATHBIO-SOLEXA2_30TUEAAXX:3:29:748:1189_length=53#0/2,SRR201779.1919192_PATHBIO-SOLEXA2_30TUEAAXX:3:29:748:1189_length=53#0/2,SRR201779.5237916_PATHBIO-SOLEXA2_30TUEAAXX:3:81:513:764_length=53#0/1,SRR201779.5237916_PATHBIO-SOLEXA2_30TUEAAXX:3:81:513:764_length=53#0/1,

# 2) normalJunctions_candidates.txt
# chr3_42701320_+:chr3_42702977_+	5	6	4	2	42701276	42703014	GT	AG	1	1	1657	100	100 	0	0	ENSG00000114853.9	ENSG00000114853.9	ZBTB47	ZBTB47	protein_coding	protein_coding	SRR201779.2086482_PATHBIO-SOLEXA2_30TUEAAXX:3:31:970:1103_length=53#0/1,SRR201779.5570454_PATHBIO-SOLEXA2_30TUEAAXX:3:86:423:861_length=53#0/2,SRR201779.5570454_PATHBIO-SOLEXA2_30TUEAAXX:3:86:423:861_length=53#0/2,SRR201779.5058273_PATHBIO-SOLEXA2_30TUEAAXX:3:79:152:958_length=53#0/2,SRR201779.652134_PATHBIO-SOLEXA2_30TUEAAXX:3:10:1049:1322_length=53#0/2,SRR201779.1924659_PATHBIO-SOLEXA2_30TUEAAXX:3:29:1578:386_length=53#0/2,


# will exit if there is an error or in a pipe
set -e -o pipefail

# In case the user does not provide any input file
###################################################
if [[ ! -e "$1" ]] || [[ ! -e "$2" ]] || [[ ! -e "$3" ]]
then
    echo "" >&2
    echo "*** ChimSplice V0.9.1 ***" >&2
    echo "" >&2
    echo "Usage:    bash ChimSplice.sh alignment_files_paths.txt genome_index.gem annot.gtf readDirectionality consensusSS outDir" >&2
    echo "" >&2
    echo "Example:  ChimSplice.sh alignment_files_paths.txt Homo_sapiens.GRCh37.chromosomes.chr.gem gencode.v7.annotation_exons.gtf UNSTRANDED '(GT,AG)' OutputDir" >&2
    echo "" >&2
    echo "Make chimeric and normal junctions from spliced aligned reads, annotated them and produce a file for each type of junction" >&2
    echo "with the number of reads supporting them, the overlapping genes and other useful info associated" >&2
    echo "" >&2
    echo "Input:" >&2
    echo "1) file listing the absolute paths to a set of BAM or GEM alignment files (GEM produced with the rna-mapper)(one row for each path). Mandatory" >&2
    echo "2) GEM indexed genome. Mandatory" >&2
    echo "3) reference annotation. Mandatory" >&2
    echo "4) readDirectionality (MATE1_SENSE|MATE2_SENSE|UNSTRANDED). Default: UNSTRANDED"  >&2
    echo "5) consensus splice sites used for splice alignments. Default='GT+AG,GC+AG,ATATC+A.,GTATC+AT'"  >&2
    echo "6) Output directory. Default: current working directory"  >&2
    echo "" >&2
    echo "Output:" >&2 
    echo "1) chimJunctions_candidates.txt" >&2
    echo "2) normalJunctions_candidates.txt" >&2
    echo "" >&2
    exit -1
fi

# In case the user does not provide any indexed genome, annotation 
##################################################################
# file or output dir or strandedness we provide default values
##############################################################
input=$1;
genome=$2;
annot=$3;

if [ ! -n "$4" ]
then
	readDirectionality="UNSTRANDED"
	spliceSites='(GT,AG),(GC,AG),(ATATC,A.),(GTATC,AT)'
	outDir=.
else
	readDirectionality=$4
	
	if [ ! -n "$5" ]
	then
		spliceSites='(GT,AG),(GC,AG),(ATATC,A.),(GTATC,AT)'
		outDir=.
	else
		spliceSites=$5
		
		if [ ! -n "$6" ]
		then
			outDir=.
		else
			outDir=$6
		fi
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

## Set bin and bash and awk scripts directories
awkDir=$rootDir/../awk
bashDir=$rootDir

# Programs and scripts
########################

## Bash
MAKE_SPLICEJUNC=$bashDir/make_spliceJunctions.sh

## Awk
SPLICEJUNC2GFF=$awkDir/spliceJunction2gff.awk
GFF2GFF=$awkDir/gff2gff.awk
JUNC_OVERLAP_DIST=$awkDir/compute_juncPercOverlap2Exon_juncDist2ExonSS.awk
SELECT_OVERLAPPING_EXONS=$awkDir/select_overlappingExons_spliceSide.awk
ADD_ANNOT_INFO=$awkDir/addAnnotInfo2spliceJunctions.awk
CLASSIFY_SPLICEJUNC=$awkDir/classifySpliceJunctions.awk


#################################################
# 1) MAKE SPLICE JUNCTIONS WITH INFO ASSOCIATED #
#################################################
# - $outDir/spliceJunc_nbStag_nbtotal_NbUnique_nbMulti_sameChrStr_okGxOrder_dist_ss1_ss2_readIds.txt

echo "1. Make splice junctions associated with their number of supporting reads, splice sites.. from split-mapping blocks" >&2

bash -f $MAKE_SPLICEJUNC $input $genome $readDirectionality $spliceSites $outDir 

################################
# 2) ANNOTATE SPLICE JUNCTIONS #
################################

echo "2. Annotate splice junctions" >&2

# 2.1) Make gff containing both parts of the junction
#####################################################
# - $outDir/spliceJunctions_2parts_beg_end.gff.gz

echo "2.1. Make a gff where each junction side correspond to an entry" >&2

awk -v OFS='\t' -f $SPLICEJUNC2GFF $outDir/spliceJunc_nbStag_nbtotal_NbUnique_nbMulti_donor_acceptor_beg_end_sameChrStr_okGxOrder_dist_readIds.txt | awk -f $GFF2GFF | gzip > $outDir/spliceJunctions_2parts_beg_end.gff.gz

# 2.2) Intersect with the annotated exons
##########################################
# - $outDir/spliceJunctions_2parts_beg_end_intersected.txt.gz

echo "2.2. Intersect each junction side with the annotated exons using bedtools" >&2

bedtools intersect -wao -a $outDir/spliceJunctions_2parts_beg_end.gff.gz -b $annot | gzip > $outDir/spliceJunctions_2parts_beg_end_intersected.txt.gz

# 2.3) For each intersection add the percentage of overlap between 
###################################################################
# block/exon and the distance to the exon boundary 
####################################################
# - $outDir/spliceJunctions_2parts_beg_end_intersected_percOverlap_distExonSS.txt.gz

echo "2.3. For each intersection add the percentage of overlap between block/exon and the distance between the junction site and the exon boundary" >&2

awk -f $JUNC_OVERLAP_DIST <( gzip -dc $outDir/spliceJunctions_2parts_beg_end_intersected.txt.gz) | gzip > $outDir/spliceJunctions_2parts_beg_end_intersected_percOverlap_distExonSS.txt.gz
 
# 2.4) Merge intersections
###########################

echo "2.4. For each splice junction side select those overlaping exons which maximize the percentage of overlap and minimize the distance to the exon boundary" >&2

awk -f $SELECT_OVERLAPPING_EXONS <( gzip -dc $outDir/spliceJunctions_2parts_beg_end_intersected_percOverlap_distExonSS.txt.gz) > $outDir/spliceJunctions_2parts_beg_end_intersected_percOverlap_distExonSS_merged.txt

# 2.5) Add annotation info to the splice junctions file with 
###############################################################

echo "2.5. Add annotation info to the splice junctions file generated in 2" >&2

awk -v juncAnnotated=$outDir/spliceJunctions_2parts_beg_end_intersected_percOverlap_distExonSS_merged.txt -f $ADD_ANNOT_INFO $outDir/spliceJunc_nbStag_nbtotal_NbUnique_nbMulti_donor_acceptor_beg_end_sameChrStr_okGxOrder_dist_readIds.txt > $outDir/spliceJunctions_annotated.txt


##################################################################################################
# 3) CLASSIFY SPLICE JUNCTIONS INTO: 															 #
#     A) NORMAL:  BOTH MATES MAPPING IN ONLY ONE GENE										     #
#     B) CHIMERIC:  BOTH MATES MAPPING IN DIFFERENT GENES										 #
#     C) UNANNOTATED: BOTH MATES MAP AND AT LEAST ONE OF THEM DO NOT OVERLAP ANY ANNOTATED GENE  #
##################################################################################################
# - $outDir/spliceJunctions_classified.txt
# - $outDir/normal_spliceJunctions.txt
# - $outDir/chimeric_spliceJunctions.txt
# - $outDir/unannotated_spliceJunctions.txt

echo "3. Classify splice junctions in normal, chimeric and unannotated" >&2

awk -v OFS='\t' -f $CLASSIFY_SPLICEJUNC $outDir/spliceJunctions_annotated.txt > $outDir/spliceJunctions_classified.txt

### Make a different file for each type of splice junction:

## A) NORMAL
awk '($2=="NORMAL")||(NR==1)' $outDir/spliceJunctions_classified.txt > $outDir/normal_spliceJunctions.txt

## B) CHIMERIC
awk '($2=="CHIMERIC")||(NR==1)' $outDir/spliceJunctions_classified.txt > $outDir/chimeric_spliceJunctions.txt

## C) UNANNOTATED
awk '($2=="UNANNOTATED")||(NR==1)' $outDir/spliceJunctions_classified.txt > $outDir/unannotated_spliceJunctions.txt


######################
# 4) CLEANUP AND END #
######################
echo "4. Cleanup and end" >&2
rm $outDir/spliceJunctions_2parts_beg_end.gff.gz $outDir/spliceJunc_nbStag_nbtotal_NbUnique_nbMulti_donor_acceptor_beg_end_sameChrStr_okGxOrder_dist_readIds.txt $outDir/spliceJunctions_2parts_beg_end_intersected.txt.gz $outDir/spliceJunctions_2parts_beg_end_intersected_percOverlap_distExonSS.txt.gz $outDir/spliceJunctions_2parts_beg_end_intersected_percOverlap_distExonSS_merged.txt $outDir/spliceJunctions_annotated.txt $outDir/spliceJunctions_classified.txt





