#!/bin/bash
# 20/03/2014

# It takes as input a chimeric junctions file and it finds the splice sites GT-AG or GC-AG and their complementary reverse AC-CT or GC-CT, adds this information to the chimeric junctions matrix and, for unstranded data, writes the junctions in biological order (first part donor and second part of the junction acceptor) and add strand information to the junction whenever it is possible. 

## Input (14 fields chimeric junctions file)
# chrX_39164738:chr10_101180417 1 2 39164722 101180450 0 NA NA ENSG00000235304.1, ENSG00000120053.8, RP11-265P11.2, GOT1, lincRNA, protein_coding,
# chr14_91739676:chr2_97494184 1 3 91739692 97494217 0 NA NA ENSG00000015133.13, ENSG00000168763.10, CCDC88C, CNNM3, protein_coding, protein_coding,

## Output (16 fields chimeric junctions file)
# chrX_39164738:chr10_101180417 1 2 39164722 101180450 0 NA NA AC CT ENSG00000235304.1, ENSG00000120053.8, RP11-265P11.2, GOT1, lincRNA, protein_coding,
# chr14_91739676:chr2_97494184 1 3 91739692 97494217 0 NA NA NA CT ENSG00000015133.13, ENSG00000168763.10, CCDC88C, CNNM3, protein_coding, protein_coding,

# usage
#######
# find_consensus_splice_sites2.sh junctions.txt strandedness outputDir

# where "stranded" is 0 to run in unstranded mode and !=0 to run it in stranded mode. If no strandedness or output directory is specified it will run in unstranded mode and the current working directory will be the output directory

# Notes
#######
# - Made for using on a 64 bit linux architecture
# - uses awk scripts


# In case the user does not provide any input file, an error message is raised
##############################################################################
if [ ! -n "$1" ]
then
    echo "" >&2
    echo Usage:    find_consensus_splice_sites2.sh junctions.txt strandedness outputDir >&2
    echo "" >&2
    echo Example:  find_consensus_splice_sites2.sh junctions.txt 0 /users/rg/brodriguez/Projects/Chimeras/Results  >&2
    echo "" >&2
    echo Takes as input a chimeric junctions file and it finds the splice sites GT-AG or GC-AG and their complementary>&2
    echo "reverse AC-CT or GC-CT, adds this information to the chimeric junctions matrix and, for unstranded data," >&2
    echo "writes the junctions in biological order (first part donor and second part of the junction acceptor)" >&2
    echo "and add strand information to the junction whenever it is possible. If no strandedness or output directory" >&2 
    echo "is specified it will run in unstranded mode and the current working directory will be the output directory" >&2 
    echo "" >&2
    exit 0
fi

input=$1

# In case the user does not provide any strandedness 
####################################################
# or output directory, default values are provided
##################################################
if [ ! -n "$2" ]
then
outdir=.
stranded=0
else
stranded=$2
if [ ! -n "$3" ]
then
outdir=.
else
outdir=$3
fi
fi

# Directories
#############
# IMPORTANT! rootDir is an environmental variable defined and exported in the main script "ChimPipe.sh" which contains the path to the root folder of ChimPipe pipeline. 
binDir=$rootDir/bin
awkDir=$rootDir/src/awk

# Programs
###########
SUBSEQ=$binDir/miscellaneous/chr_subseq 
JUNCBIOL=$awkDir/junctions_bio_order.awk
GFF2GFF=$awkDir/gff2gff.awk

# Setting input variables
#########################
input=$1
stranded=$2
outdir=$3

# Subtract the two pairs of nucleotides that should correspond to the donor and acceptor respectively. 
######################################################################################################
# If they are not the cannonical writes NA
##########################################
awk '{print $1}' $input | while read junc  
do 
	echo $junc | awk '{split($1,a,":");split(a[1],b1,"_");split(a[2],b2,"_"); print b1[1], b1[2], b2[1], b2[2]}' | while read chr1 coord1 chr2 coord2 
	do
		dnt1=`$SUBSEQ /db/seq/genomes/H.sapiens/golden_path_200902/chromFa/$chr1.fa $(($coord1 + 1)) $(($coord1 + 2)) | sed 's/.*/\U&/' | awk '{if (($1=="GT")||($1=="AG")||($1=="AC")||($1=="CT")||($1=="GC")){print $1;}else{print "NA";}}'`
		dnt2=`$SUBSEQ /db/seq/genomes/H.sapiens/golden_path_200902/chromFa/$chr2.fa $(($coord2 - 2)) $(($coord2 - 1)) | sed 's/.*/\U&/' | awk '{if (($1=="GT")||($1=="AG")||($1=="AC")||($1=="CT")||($1=="GC")){print $1;}else{print "NA";}}'`
		echo $junc $dnt1 $dnt2 
	done
done > $outdir/junct_dnt1_dnt2.txt


# Adds the splice sites information to the matrix and for unstranded data writes the junctions in biological order 
##################################################################################################################
# and add strand information to the junction whenever it is possible. 
#####################################################################
awk -v fileRef=$outdir/junct_dnt1_dnt2.txt -v stranded=$stranded -f $JUNCBIOL $input 

# Cleaning. 
###########
rm $outdir/junct_dnt1_dnt2.txt


