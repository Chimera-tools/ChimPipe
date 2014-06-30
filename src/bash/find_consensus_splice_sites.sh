#!/bin/bash


# It takes as input a chimeric junctions file and the indexed genome and it finds the splice sites GT-AG or GC-AG and their complementary reverse AC-CT or GC-CT, adds this information to the chimeric junctions matrix and, for unstranded data, writes the junctions in biological order (first part donor and second part of the junction acceptor) and add strand information to the junction whenever it is possible. 

## Input (14 fields chimeric junctions file)
# chrX_39164738:chr10_101180417 1 2 39164722 101180450 0 NA NA ENSG00000235304.1, ENSG00000120053.8, RP11-265P11.2, GOT1, lincRNA, protein_coding,
# chr14_91739676:chr2_97494184 1 3 91739692 97494217 0 NA NA ENSG00000015133.13, ENSG00000168763.10, CCDC88C, CNNM3, protein_coding, protein_coding,

## Output (16 fields chimeric junctions file)
# chrX_39164738:chr10_101180417 1 2 39164722 101180450 0 NA NA AC CT ENSG00000235304.1, ENSG00000120053.8, RP11-265P11.2, GOT1, lincRNA, protein_coding,
# chr14_91739676:chr2_97494184 1 3 91739692 97494217 0 NA NA NA CT ENSG00000015133.13, ENSG00000168763.10, CCDC88C, CNNM3, protein_coding, protein_coding,

# usage
#######
# find_consensus_splice_sites2.sh junctions.txt indexed_genome strandedness outputDir

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
    echo Usage:    find_consensus_splice_sites2.sh junctions.txt genome.gem strandedness outputDir >&2
    echo "" >&2
    echo Example:  find_consensus_splice_sites2.sh junctions.txt /users/rg/sdjebali/ENCODE_AWG/Analyses/Mouse_Human/Chimeras/Human/Homo_sapiens.GRCh37.chromosomes.chr.gem 0 /users/rg/brodriguez/Projects/Chimeras/Results  >&2
    echo "" >&2
    echo "Takes as input a chimeric junctions file, an indexed genome, strandedness (1 stranded and 0 unstranded) and it finds the">&2
    echo "splice sites GT-AG or GC-AG and their complementary reverse AC-CT or GC-CT, adds this information to the chimeric junctions " >&2
    echo "matrix and, for unstranded data, writes the junctions in biological order (first part donor and second part of the junction" >&2
    echo " acceptor) and add strand information to the junction whenever it is possible. If no strandedness or output directory" >&2 
    echo "is specified it will run in unstranded mode and the current working directory will be the output directory" >&2 
    echo "" >&2
    exit 0
fi
input=$1
genome=$2

# In case the user does not provide any strandedness 
####################################################
# or output directory, default values are provided
##################################################
if [ ! -n "$3" ]
then
outdir=.
stranded=0
else
stranded=$3
if [ ! -n "$4" ]
then
outdir=.
else
outdir=$4
fi
fi

# Directories 
#############
# Environmental variables 
# rootDir - path to the root folder of ChimPipe pipeline. 
# It is an environmental variable defined and exported in the main script 
binDir=$rootDir/bin
awkDir=$rootDir/src/awk

# Programs
###########
SUBSEQ=$binDir/miscellaneous/chr_subseq 
JUNCBIOL=$awkDir/junctions_bio_order.awk
GFF2GFF=$awkDir/gff2gff.awk
RETRIEVER=$binDir/gemtools-1.7.1-i3/bin/gem-retriever

# Subtract the two pairs of nucleotides that should correspond to the donor and acceptor respectively. 
######################################################################################################
# If they are not the cannonical writes NA
##########################################

awk '{print $1}' $input | while read junc  
do 
	echo $junc | awk -v OFS='\t' '{split($1,a,":");split(a[1],b1,"_");split(a[2],b2,"_");print b1[1], "+", (b1[2]+1), (b1[2]+2); print b2[1], "+", (b2[2]-2), (b2[2]-1)}'
done | $RETRIEVER $genome | sed 's/.*/\U&/' |  awk '{if (($1=="GT")||($1=="AG")||($1=="AC")||($1=="CT")||($1=="GC")){print $1;}else{print "NA";}}' > $outdir/consensusSS.txt

awk -v consensusSS=$outdir/consensusSS.txt 'BEGIN{junc=1;part=1;while(getline < consensusSS >0){ss[junc]=ss[junc]" "$1; if (part==2){part=1;junc++}else{part++}}junc=0}{junc++; split(ss[junc],a," ");print $1, a[1], a[2]}' $input > $outdir/junct_ss1_ss2.txt


# Adds the splice sites information to the matrix and for unstranded data writes the junctions in biological order 
##################################################################################################################
# and add strand information to the junction whenever it is possible. 
#####################################################################
awk -v fileRef=$outdir/junct_ss1_ss2.txt -v stranded=$stranded -f $JUNCBIOL $input   

# Cleaning. 
###########
rm $outdir/consensusSS.txt $outdir/junct_ss1_ss2.txt


