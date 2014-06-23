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
awk -v fileRef=$outdir/junct_dnt1_dnt2.txt -v stranded=$stranded 'BEGIN{while (getline < fileRef >0){part1[$1]=$2;part2[$1]=$3;if (($2=="GT")||($2=="AG")){strand1[$1]="+";}else{if (($2=="AC")||($2=="CT")){strand1[$1]="-";}else{strand1[$1]=".";}}if (($3=="GT")||($3=="AG")){strand2[$1]="+";}else{if (($3=="AC")||($3=="CT")){strand2[$1]="-";}else{strand2[$1]=".";}}}}{split($1,a,":");split(a[1],a1,"_");split(a[2],a2,"_");if ((((part1[$1]=="AG")&&((part2[$1]=="GT")||(part2[$1]=="AC")||(part2[$1]=="GC"))) || ((part1[$1]=="CT")&&((part2[$1]=="GT")||(part2[$1]=="AC")||(part2[$1]=="GC")))) && (stranded==0)){junc=a2[1]"_"a2[2]"_"strand2[$1]":"a1[1]"_"a1[2]"_"strand1[$1];ss1=part2[$1];ss2=part1[$1];gnlist1=$10;gnlist2=$9;gnname1=$12;gnname2=$11;bt1=$14;bt2=$13;}else{if (stranded==0){junc=a1[1]"_"a1[2]"_"strand1[$1]":"a2[1]"_"a2[2]"_"strand2[$1];}else{junc=a1[1]"_"a1[2]"_"a1[3]":"a2[1]"_"a2[2]"_"a2[3];}ss1=part1[$1];ss2=part2[$1];gnlist1=$9;gnlist2=$10;gnname1=$11;gnname2=$12;bt1=$13;bt2=$14;}print junc, $2, $3, $4, $5, $6, $7, $8, ss1, ss2, gnlist1, gnlist2, gnname1, gnname2, bt1, bt2}' $input > $outdir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt 

# Cleaning. 
###########
rm $outdir/junct_dnt1_dnt2.txt


