#!/bin/bash

<<authors
*****************************************************************************
	
	make_spliceJunctions.sh
	
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

# Takes as input a gzipped GFF file with splice split-mapping blocks and makes the correct splice junctions joining the donnor and acceptor sides of the blocks, compute the number of staggered reads, the total number of reads, the number of unique and multi mappings supporting the junction, their maximum beginning and end coordinates, their consensus splice sites.... It does it in three consecutive steps:

# 1) Make staggered reads from split-mapping blocks
# 2) Substract with the gem retriever 8-kmer on each side of the blocks. This will produce 4 8-kmer, 2 of them containing the splice site sequence. 
# 3) Make the chimeric junctions looking at the splice sites in the k-mers. 

## Input (split-alignment blocks)
# chr19	ChimPipe	alBlock1	54677696	54677731	1	-	.	ReadName: "SRR201779.884436_PATHBIO-SOLEXA2_30TUEAAXX:3:13:198:2046_length=53#0/2";
# chr3	ChimPipe	alBlock2	41547778	41547794	1	+	.	ReadName: "SRR201779.884436_PATHBIO-SOLEXA2_30TUEAAXX:3:13:198:2046_length=53#0/2";
# chr10	ChimPipe	alBlock1	124542237	124542251	1	+	.	ReadName: "SRR201779.884531_PATHBIO-SOLEXA2_30TUEAAXX:3:13:1048:1883_length=53#0/1";
# chrX	ChimPipe	alBlock2	66765937	66765974	1	+	.	ReadName: "SRR201779.884531_PATHBIO-SOLEXA2_30TUEAAXX:3:13:1048:1883_length=53#0/1";

## Output (splice junctions)
# chr11_60781486_+:chr11_120100337_+ 1 1 60781475 120100375 GT AG
# chr16_30537851_-:chr15_40629945_+ 1 1 30537885 40629960 AC AG
# chr13_64414915_-:chr18_74690910_- 1 3 64414930 74690876 AC CT
# chr11_71714505_+:chr22_50928336_- 1 1 71714472 50928320 GT CT

# usage
#######
# make_chimjunctions.sh blocks.gff genome.gem strandedness spliceSites outputDir

# where "stranded" is 0 to run in unstranded mode and !=0 to run it in stranded mode. If no strandedness or output directory is specified it will run in unstranded mode and the current working directory will be the output directory

# Notes
#######
# - Made for using on a 64 bit linux architecture
# - uses awk scripts

# will exit if there is an error or in a pipe
set -e -o pipefail

# In case the user does not provide any input file, an error message is raised
##############################################################################
if [ ! -n "$1" ] || [ ! -n "$2" ]
then
    echo "" >&2
    echo Usage:  make_spliceJunctions.sh input_alignments.txt genome.gem strandedness spliceSites outputDir >&2
    echo "" >&2
    echo Example:  make_spliceJunctions.sh spliceAlignments_2blocks.gff.gz Homo_sapiens.GRCh37.chromosomes.chr.gem 0 'GT+AG' ./Results  >&2
    echo "" >&2
    exit -1
fi
input=$1
genome=$2

# In case the user does not provide any strandedness 
####################################################
# or output directory, default values are provided
##################################################
if [ ! -n "$3" ]
then
	spliceSites='GT+AG,GC+AG,ATATC+A.,GTATC+AT'
	outDir=.
	readDirectionality="UNSTRANDED"
else
	readDirectionality=$3
	if [ ! -n "$4" ]
	then
		spliceSites='GT+AG,GC+AG,ATATC+A.,GTATC+AT'
		outDir=.
	else
		spliceSites=$4
		if [ ! -n "$5" ]
		then
			outDir=.
		else
			outDir=$5
		fi
	fi
fi

if [[ "$readDirectionality" == "MATE1_SENSE" ]] || [[ "$readDirectionality" == "MATE2_SENSE" ]]
then
	stranded=1;
elif [ "$readDirectionality" == "UNSTRANDED" ]
then
	stranded=0;
else 
	exit;
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

## Set bin and awk directories
binDir=$rootDir/../../bin
awkDir=$rootDir/../awk

# Programs
###########
## Awk
BED2BEDPE=$awkDir/bed2bedPE.awk
BEDPE_CORRECTSTRAND=$awkDir/bedPECorrectStrand.awk
BEDPE2GFF=$awkDir/bedPE2gff.awk
GFF2GFF=$awkDir/gff2gff.awk
GEM_CORRECT_STRAND=$awkDir/gemCorrectStrand.awk
GEM2GFF=$awkDir/gemsplit2gff.awk
MAKE_STAGGERED=$awkDir/make_staggeredSplitMappings.awk

MAKE_JUNCTIONS=$awkDir/staggered2junct.awk

## Bin
RETRIEVER=$binDir/gemtools-1.7.1-i3/gem-retriever

# 1. Make a gff with splice junction blocks from BAM and GEM alignment files
################################################################################
# - $outDir/spliceAlignments_2blocks.gff.gz

nbFiles=1; # Count files

### Iterate reading a BAM or GEM file in each iteration and saving junction blocks into a GFF
cat $input | while read file
do
	fileName=$(basename "$file")
	extension="${fileName##*.}"
	
	## A) First file, create output gff file
	if [[ "$nbFiles" == "1" ]];
    then
		case $extension in
    	"bam")
       		samtools view -bF 4 $file | bedtools bamtobed -bed12 -tag NH -i stdin | awk '($10==2) && ($1 !~ /M/) && ($1 !~ /Mt/) && ($1 !~ /MT/)' | awk -v rev='1' -f $BED2BEDPE | awk -v readDirectionality=$readDirectionality -f $BEDPE_CORRECTSTRAND | awk -f $BEDPE2GFF | awk -f $GFF2GFF  > $outDir/spliceAlignments_2blocks.gff
        	;;	    
    	"map")
    	    awk -v readDirectionality=$readDirectionality -f $GEM_CORRECT_STRAND $file | awk '($5 !~ /M/) && ($5 !~ /Mt/) && ($5 !~ /MT/)' | awk -v rev="0" -f $GEM2GFF | awk -f $GFF2GFF > $outDir/spliceAlignments_2blocks.gff
            ;;    
   		 esac
   	
   	## B) Not first file, append to already created output gff file
	else
		case $extension in
    	"bam")
       		samtools view -bF 4 $file | bedtools bamtobed -bed12 -tag NH -i stdin | awk '($10==2) && ($1 !~ /M/) && ($1 !~ /Mt/) && ($1 !~ /MT/)' | awk -v rev='1' -f $BED2BEDPE | awk -v readDirectionality=$readDirectionality -f $BEDPE_CORRECTSTRAND | awk -f $BEDPE2GFF | awk -f $GFF2GFF  >> $outDir/spliceAlignments_2blocks.gff
        	;;	    
    	"map")
    	    awk -v readDirectionality=$readDirectionality -f $GEM_CORRECT_STRAND $file | awk '($5 !~ /M/) && ($5 !~ /Mt/) && ($5 !~ /MT/)' | awk -v rev="0" -f $GEM2GFF | awk -f $GFF2GFF >> $outDir/spliceAlignments_2blocks.gff
            ;;
         esac
	fi
	
	# Update counter 
	let nbFiles=nbFiles+1;
done 

### Compress gff file
gzip $outDir/spliceAlignments_2blocks.gff


# 2. Make the staggered reads from gff containing split-mapping blocks
######################################################################## 
zcat $outDir/spliceAlignments_2blocks.gff.gz | awk -f $MAKE_STAGGERED > $outDir/staggered_nbTotal_nbUnique_nbMulti_readIds.txt


# 3. Subtract for each staggered the first 8 flanking nucleotides in each side (this makes a total of 4 8-mers). 
###############################################################################################################
paste <(cat $outDir/staggered_nbTotal_nbUnique_nbMulti_readIds.txt) <(awk -v OFS='\t' '{split($1,a,":"); split(a[1],a1,"_"); chrA=a1[1]; begA=a1[2]; endA=a1[3]; strA=a1[4]; split(a[2],a2,"_"); chrB=a2[1]; begB=a2[2]; endB=a2[3]; strB=a2[4]; print chrA, strA, (begA-8), (begA-1); print chrA, strA, (endA+1), (endA+8); print chrB, strB, (begB-8), (begB-1); print chrB, strB, (endB+1), (endB+8)}' $outDir/staggered_nbTotal_nbUnique_nbMulti_readIds.txt | $RETRIEVER $genome | sed 's/.*/\U&/' | awk 'BEGIN{counter=1;key=1;}{seq[key]=seq[key]" "$1; counter++; if (counter==5) {counter=1; key++}}END{for (i=1;i<key;i++){print seq[i]}}') | awk '{print $1, $2, $3, $4, $6, $7, $8, $9, $5;}' | awk -f $GFF2GFF > $outDir/staggered_nbTotal_nbUnique_nbMulti_nts_readIds.txt

# 4. Make the splice junctions, compute the number of staggered and total reads supporting the junctions,  
########################################################################################################
# the consensus splice sites,  the maximum beginning and end plus some further info..
#######################################################################################
awk -v strandedness=$stranded -v CSS=$spliceSites -f $MAKE_JUNCTIONS $outDir/staggered_nbTotal_nbUnique_nbMulti_nts_readIds.txt > $outDir/spliceJunc_nbStag_nbtotal_NbUnique_nbMulti_sameChrStr_okGxOrder_dist_ss1_ss2_readIds.txt


# Cleaning. 
###########
rm $outDir/spliceAlignments_2blocks.gff.gz $outDir/staggered_nbTotal_nbUnique_nbMulti_readIds.txt $outDir/staggered_nbTotal_nbUnique_nbMulti_nts_readIds.txt 



