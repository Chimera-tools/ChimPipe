#!/bin/bash
# It takes as input a chimeric junctions file and it finds the splice sites GT-AG or GC-AG and their complementary reverse AC-CT or GC-CT. It adds this information to the matrix and for unstranded data writes the junctions in biological order (first part donor and second part of the junction acceptor) and add strand information to the junction whenever it is possible. 

# Input (14 fields chimeric junctions file)
# chrX_39164738:chr10_101180417 1 2 39164722 101180450 0 NA NA ENSG00000235304.1, ENSG00000120053.8, RP11-265P11.2, GOT1, lincRNA, protein_coding,
# chr14_91739676:chr2_97494184 1 3 91739692 97494217 0 NA NA ENSG00000015133.13, ENSG00000168763.10, CCDC88C, CNNM3, protein_coding, protein_coding,

# usage
#######
# find_consensus_splice_sites.sh junctions.txt stranded outdir

# where "stranded" is 0 to run in unstranded mode and !=0 to run it in stranded mode

# Notes
#######
# - Made for using on a 64 bit linux architecture
# - uses awk scripts

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
awk -v fileRef=$outdir/junct_dnt1_dnt2.txt -v stranded=$stranded -f $JUNCBIOL $input > $outdir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt 

# Cleaning. 
###########
rm $outdir/junct_dnt1_dnt2.txt


