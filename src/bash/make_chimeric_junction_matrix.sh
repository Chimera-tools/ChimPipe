#!/bin/bash

<<authors
*****************************************************************************
	
	make_chimeric_junction_matrix.sh
	
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

# usage
#######
# make_chimeric_junction_matrix.sh junctions annot stranded outdir

# where:
#		junctions.  Input file containing the junctions
# 		annot. 		Annotation file in gtf format 
# 		stranded. 	Data strandedness. 0 for unstranded and 1 for stranded data
# 		Outdir. 	Output directory. 

# Input (chimeric junctions file):
# chrX_39164738_+:chr10_101180417_+ 1 2 39164722 101180450
# chr14_91739676_.:chr2_97494184_- 1 3 91739692 97494217
# chr16_69752040_+:chr14_101014151_+ 1 1 69752026 101014116

# Output: 
# chrX_39164738_+:chr10_101180417_+ 1 2 39164722 101180450 0 NA NA ENSG00000235304.1, ENSG00000120053.9, RP11-265P11.2, GOT1, lincRNA, protein_coding,
# chr14_91739676_.:chr2_97494184_- 1 3 91739692 97494217 0 NA NA ENSG00000015133.14, ENSG00000168763.11, CCDC88C, CNNM3, protein_coding, protein_coding,
# chr16_69752040_+:chr14_101014151_+ 1 1 69752026 101014116 0 NA NA ENSG00000181019.8, ENSG00000183092.11, NQO1, BEGAIN, protein_coding, protein_coding,

# Notes
#######
# - Made for using on a 64 bit linux architecture
# - uses awk scripts

# Directories 
#############
# Environmental variables 
# rootDir - path to the root folder of ChimPipe pipeline. 
# It is an environmental variable defined and exported in the main script 
binDir=$rootDir/bin
awkDir=$rootDir/src/awk

# Programs
###########
CUTGFF=$awkDir/cutgff.awk
GFF2GFF=$awkDir/gff2gff.awk
JUNC2GFF=$awkDir/chimjunc2gff.awk 
OVERLAP=$binDir/overlap

# Setting input variables
#########################
input=$1
annot=$2
stranded=$3
outdir=$4

# In case the user does not provide any input file
###################################################
if [ ! -n "$1" ]
then
    echo "" >&2
    echo Usage: make_chimeric_junction_matrix2.sh file_with_junc_one_sample.txt annot outdir>&2
    echo "" >&2
    echo Example: make_chimeric_junction_matrix.sh ~brodriguez/Chimeras_project/Chimeras_detection_pipeline/Chimera_mapping/Benchmark/Results/V0.5.0/Negative_set_calogero/lib50_Q2/UCSCgenesCalogero/Chimsplice/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt /users/rg/projects/encode/scaling_up/whole_genome/Gencode/version10/Long/Element/gen10.long.exon.gtf ~brodriguez/Analysis >&2
    echo "" >&2
    echo Takes as input:
    echo 1\) a matrix with with chimeric junctions, nb stag, nb total, max beg & max end >&2
    echo 2\) a gff annotation file with at least exon rows associated to their gene id in column 10 >&2
    echo 3\) Data strandedness. 0 for unstranded and 1 for stranded data. If not defined it will run in unstranded mode >&2
    echo 4\) Output directory. If not defined it will be the current working directory >&2
    echo and add to this matrix information about their position, distance and overlapping genes>&2 
    echo "" >&2
    exit -1
else
input=$1
fi

# In case the user does not provide any annotation file
########################################################
# or output dir or strandedness we provide default values
#########################################################
if [ ! -n "$2" ]
then
    annot=/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version19/Long/gencode.v19.annotation.long.gtf
	stranded=0
	outdir="."
else
    annot=$2
    if [ ! -n "$3" ]
    then
	stranded=0
	outdir="."
    else
	stranded=$3
		if [ ! -n "$4" ]
    	then
		outdir="."
    	else
		outdir=$4
    	fi
    fi
fi

# Display script configuration
##############################

#printf "\n"
#printf "*****Make chimeric junction matrix configuration*****\n"
#printf "Input: $input\n"
#printf "Annotation: $annot\n"
#printf "Strand: $stranded\n"
#printf "Outdir: $outdir\n"
#printf "\n"
#exit 1



# 0) Make an exon file with 20 fields formatted ok for overlap and determine the fields of 
##########################################################################################
# the number of the column where gene_name and gene_type information are located
################################################################################
echo I am making the exon gff file and determining where the gene_name and gene_type lie \in this gff >&2
annotbase=`basename ${annot%.gtf}`
awk '$3=="exon"' $annot | awk -v to=20 -f $CUTGFF > $outdir/$annotbase.exons.20flds.gff
str=`awk 'NR==10{fldgnname="."; fldgnbt="."; k=9; while(k<=(NF-1)){if($k=="gene_name"){fldgnname=(k+1);}else{if($k=="gene_type"){fldgnbt=(k+1);}} k+=2}}END{print fldgnname, fldgnbt}' $outdir/$annotbase.exons.20flds.gff`
fldgnname=`echo $str | awk '{print $1}'` 
fldgnbt=`echo $str | awk '{print $2}'`
echo they are $fldgnname and $fldgnbt respectively >&2

# 1) Add 2 booleans, one for intrachrsamestr and one for ok genomic order when same intrachrsamestr 
###################################################################################################
# as well as the distance when it makes sense
#############################################
echo I am adding the information of intrachrstr, okgxorder and distance when it makes sense >&2

awk '{split($1,a,":"); split(a[1],a1,"_"); split(a[2],a2,"_"); samechrstr=(((a1[1]==a2[1])&&(a1[3]==a2[3])) ? 1 : 0); okgxorder=((samechrstr==1) ? ((((a1[3]=="+")&&((a2[2]-a1[2])>=0))||((a1[3]=="-")&&((a1[2]-a2[2])>=0))) ? 1 : 0) : "NA"); print  $1, $2, $3, $4, $5, samechrstr, okgxorder, ((okgxorder==1) ? ((a1[3]=="+") ? (a2[2]-a1[2]) : (a1[2]-a2[2])) : "NA"), $6, $7;}' $input > $outdir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt

# 2) add the list of gene ids, real gene names, gene biotypes for each part of each junction
############################################################################################
echo I am adding the list of gene ids, gene names and gene biotypes for each part of the junction >&2
# a. first make the gff file of the chimeric junctions
######################################################
echo I am making the gff file of the junctions >&2
awk -v stranded=$stranded -f $JUNC2GFF $outdir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt | awk -f $GFF2GFF > $outdir/distinct_junctions_withmaxbegandend_two_parts.gff

# b. then overlap to know which genes have their exons overlapping each part of the junction 
#############################################################################################
echo I am overlapping the two parts of each junction with the annotated exons>&2
$OVERLAP $outdir/distinct_junctions_withmaxbegandend_two_parts.gff $outdir/$annotbase.exons.20flds.gff -st $stranded -m 10 -f ex -nr -o $outdir/distinct_junctions_withmaxbegandend_two_parts_withgnlist.gff 

# c. then add the gene real names and the gene biotypes 
#######################################################
echo I am adding the gene real names and biotypes >&2

awk -v gnname=$fldgnname -v gnbt=$fldgnbt -v fileRef=$outdir/$annotbase.exons.20flds.gff 'BEGIN{while (getline < fileRef >0){if(gnbt!="."){split($gnbt,a,"\""); bt[$10]=a[2]}else{bt[$10]="."} if(gnname!="."){split($gnname,b,"\""); name[$10]=b[2]}else{name[$10]="."}}} {s1=""; s2=""; split($NF,a,","); if(gnbt!="."){k=1; while(a[k]!=""){s1=(s1)(bt[a[k]])(",");k++}}else{s1=".";}if(gnname!="."){k=1; while(a[k]!=""){s2=(s2)(name[a[k]])(","); k++;}}else{s2=".";} print $0, "gnnamelist:", s2, "btlist:", s1}' $outdir/distinct_junctions_withmaxbegandend_two_parts_withgnlist.gff | awk -f $GFF2GFF > $outdir/distinct_junctions_withmaxbegandend_two_parts_withgnlist_gnnamelist_btlist.gff

# d. finally introduce all this information in the matrix
#########################################################
echo I am putting this information back to the matrix >&2

awk -v fileRef=$outdir/distinct_junctions_withmaxbegandend_two_parts_withgnlist_gnnamelist_btlist.gff 'BEGIN{while (getline < fileRef >0){split($10,a,":");split(a[1],a1,"_");split(a[2],a2,"_");s="";split($14,b,",");k=1;while(b[k]!=""){split(b[k],c,"\"");s=(s)(c[2])(",");k++;}if ((a1[2]==$4) || (a1[2]==$5)){gnlist[$10,1]=s;gnnamelist[$10,1]=$16;btlist[$10,1]=$18;}if ((a2[2]==$4) || (a2[2]==$5))	{gnlist[$10,2]=s;gnnamelist[$10,2]=$16;btlist[$10,2]=$18}}}{print $0, gnlist[$1,1], gnlist[$1,2], gnnamelist[$1,1], gnnamelist[$1,2], btlist[$1,1], btlist[$1,2]}' $outdir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt 

# 3) Clean
##########
echo I am cleaning >&2
rm $outdir/$annotbase.exons.20flds.gff $outdir/distinct_junctions_withmaxbegandend_two_parts.gff $outdir/distinct_junctions_withmaxbegandend_two_parts_withgnlist_gnnamelist_btlist.gff $outdir/distinct_junctions_withmaxbegandend_two_parts_withgnlist.gff $outdir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt

echo I am done >&2
