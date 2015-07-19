#!/bin/bash

<<authors
*****************************************************************************
	
	make_chimjunctions.sh
	
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

# Takes as input a file containing staggered reads split-mapping in two exons from different genes. It uses this information to make the correct chimeric junctions joining the donnor and acceptor sides of the blocks, compute the number of staggered reads supporting the junctions, their maximum beginning and end coordinates and their consensus splice sites. It does it in two consecutive steps:

# 1-) Substract with the gem retriever 8-kmer on each side of the blocks. This will will 4 8-kmer, 2 of them containing the splice site sequence. 
# 2-) Make the chimeric junctions looking at the splice sites in the k-mers. 

## Input 
# chr10_100018895_100018913_-:chr3_58279376_58279407_+ 1
# chr10_100143466_100143480_+:chr3_183901336_183901371_+ 1
# chr10_100154777_100154810_+:chr19_36214052_36214068_+ 1
# chr10_100171018_100171033_+:chr7_72419139_72419173_- 2

## Output 
# chr11_60781486_+:chr11_120100337_+ 1 1 60781475 120100375 GT AG
# chr16_30537851_-:chr15_40629945_+ 1 1 30537885 40629960 AC AG
# chr13_64414915_-:chr18_74690910_- 1 3 64414930 74690876 AC CT
# chr11_71714505_+:chr22_50928336_- 1 1 71714472 50928320 GT CT

# usage
#######
# make_chimjunctions.sh blocks.txt indexed_genome strandedness spliceSites outputDir

# where "stranded" is 0 to run in unstranded mode and !=0 to run it in stranded mode. If no strandedness or output directory is specified it will run in unstranded mode and the current working directory will be the output directory

# Notes
#######
# - Made for using on a 64 bit linux architecture
# - uses awk scripts

# In case the user does not provide any input file, an error message is raised
##############################################################################
if [ ! -n "$1" ] || [ ! -n "$2" ]
then
    echo "" >&2
    echo Usage:  make_chimjunctions.sh staggered.txt genome.gem strandedness spliceSites outputDir >&2
    echo "" >&2
    echo Example:  make_chimjunctions.sh staggered.txt /users/rg/sdjebali/ENCODE_AWG/Analyses/Mouse_Human/Chimeras/Human/Homo_sapiens.GRCh37.chromosomes.chr.gem 0 'GT+AG' /users/rg/brodriguez/Projects/Chimeras/Results  >&2
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
	outdir=.
	stranded=0
else
	stranded=$3
	if [ ! -n "$4" ]
	then
		spliceSites='GT+AG,GC+AG,ATATC+A.,GTATC+AT'
		outdir=.
	else
		spliceSites=$4
		if [ ! -n "$5" ]
		then
			outdir=.
		else
			outdir=$5
		fi
	fi
fi

# Directories 
#############
rootDir=/nfs/users/rg/brodriguez/Chimeras_project/Chimeras_detection_pipeline/ChimPipe
binDir=$rootDir/bin
awkDir=$rootDir/src/awk

# Programs
###########
RETRIEVER=$binDir/gemtools-1.7.1-i3/gem-retriever
GFF2GFF=$awkDir/gff2gff.awk
MAKE_JUNCTIONS=$awkDir/staggered2junct.awk


# Subtract the two pairs of nucleotides that should correspond to the donor and acceptor respectively. 
######################################################################################################
paste <(cat $input) <(awk -v OFS='\t' '{split($1,a,":"); split(a[1],a1,"_"); chrA=a1[1]; begA=a1[2]; endA=a1[3]; strA=a1[4]; split(a[2],a2,"_"); chrB=a2[1]; begB=a2[2]; endB=a2[3]; strB=a2[4]; print chrA, strA, (begA-8), (begA-1); print chrA, strA, (endA+1), (endA+8); print chrB, strB, (begB-8), (begB-1); print chrB, strB, (endB+1), (endB+8)}' $input  | $RETRIEVER $genome | sed 's/.*/\U&/' | awk 'BEGIN{counter=1;key=1;}{seq[key]=seq[key]" "$1; counter++; if (counter==5) {counter=1; key++}}END{for (i=1;i<key;i++){print seq[i]}}') | awk -f $GFF2GFF > $outdir/staggered_nbReads_nts.txt

# Make the chimeric junctions, compute the number of staggered and total reads supporting the junctions,  
########################################################################################################
# the consensus splice sites and the maximum beginning and end
##############################################################
awk -v strandedness=$stranded -v CSS=$spliceSites -f $MAKE_JUNCTIONS $outdir/staggered_nbReads_nts.txt

# Cleaning. 
###########
rm $outdir/staggered_nbReads_nts.txt


