#!/bin/bash

# /users/rg/brodriguez/Chimeras_project/Chimeras_detection_pipeline/Chimera_mapping/Workdir/Bash/make_chimeric_junction_matrix4.sh
##################################################################################################################################
# this script takes as input the chimeric junctions obtained with chimpipe version 0.6.6
# (ie from find_exon_exon_connections_from_splitmappings_better2.sh followed by 
# find_chimeric_junctions_from_exon_to_exon_connections_better2.sh and followed by PE information
# and nt similarity) 
# and produces a matrix with all the chimeric junctions with the strand known for the two parts (either
# guessed from splice sites for unstranded or known for strandeed) with at least 1 pe read and at least
# x stag in one exp, obtained for all exp as well as their beg and end across all experiments, 
# whether the two parts are on the same chr and strand, whether they are on the expected 
# genomic order, their distance, the list of genes which exons overlap each part of the junction, 
# the real names of those genes, the biotypes of those genes and the number of staggered reads
# supporting the junction in each exp.
# This script differs from make_chimeric_junction_matrix2.sh because it takes as input a matrix with one more
# column = the list of gene pairs with number of paired end read for each chimeric junction
# this file is not located in the chimsplice dir but one dir above

# !!! should not be run twice in the same dir at the same time (common files) !!!
# !!! the header might be wrong it would be better to take it from input file !!!

# More precisely, this script takes as input:
#############################################
# 1) a 2 column file containing the chimeric junction files associated to each exp
#   = where column 1 is a unique identifier for the exp, column 2 is the absolute
#   path to the chimeric junction file with the number of staggered reads
#   example
#   #######
#   T0 /users/rg/sdjebali/T47D/Chimeras/T0/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_PEinfo_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt
#   4 (2 fields)
#   chr3_49321897_-:chr12_113336328_+ 1 4 49321931 113336342 0 NA NA NA NA ENSG00000114316.8, ENSG00000089169.10, USP4, RPH3A, protein_coding, protein_coding, .
#   4258 (17 fields)
# 2) a minimum number of staggered reads for a junction to be reported in the matrix (e.g. 2)
# and this script will produce as output a matrix like this one where the script has been run
#############################################################################################
# junc_id beg end samechrstr okgxorder dist gnlist1 gnlist2 gnname1 gnname2 gnbt1 gnbt2 LID8963 LID8964 LID8965 LID8966 LID8969 LID8970 LID44498 LID44499 LID16629 LID16630 LID8461 LID8462 LID16633 LID16634 LID16635 LID16636 LID8710 LID8711 LID8463 LID8464 LID45016 LID45017 LID16627 LID16628 LID8686 LID8687 LID44594 LID44497 LID16631 LID16632 LID8692 LID8701 LID46598 LID46599 LID8967 LID8968
# chr17_5405057_-:chr17_5403432_- 5405127 5403366 1 1 1625 ENSG00000091592.8, ENSG00000170233.3, NLRP1, NLRP1, protein_coding, ambiguous_orf, 0 0 0 0 0 1 14 11 0 2 0 0 0 0 0 0 0 0 2 1 2 0 4 2 0 0 47 9 2 3 0 0 3 3 1 0


# usage
#######
# make_chimeric_junction_matrix3.sh file_with_junc_each_exp.txt minstag

# example
#########
# cd ~sdjebali/T47D/Chimeras/AllCases
# time make_chimeric_junction_matrix3.sh lid_junctionfile.txt 2 2> make_chimeric_junction_matrix3.err

# Notes
#######
# - Made for using on a 64 bit linux architecture
# - uses awk scripts

# Programs
###########
CUTGFF=~sdjebali/Awk/cutgff.awk
BEGEND=~sdjebali/Awk/add_maxbegend_to_chimjunc2.awk
GFF2GFF=~sdjebali/Awk/gff2gff.awk
OVERLAP=~sdjebali/bin/overlap

# In case the user does not provide any input file
###################################################
if [ ! -n "$1" ]
then
    echo "" >&2
    echo Usage: make_chimeric_junction_matrix3.sh file_with_junc_each_exp.txt minstag >&2
    echo "" >&2
    echo Example: make_chimeric_junction_matrix3.sh lid_junctionfile.txt 2 >&2
    echo "" >&2
    echo Takes as input:
    echo 1\) a 2 column file containing the chimeric junction files associated to each exp >&2
    echo where column 1 is a unique identifier for the exp, and column 2 is the absolute >&2
    echo path to the chimeric junction file with the number of staggered reads >&2
    echo 2\) a minimum number of staggered reads for a junction to be reported in the matrix >&2
    echo and produces a matrix containing the chimeric junctions obtained across all experiments >&2
    echo with additional information about their position, distance, genes and expression >&2
    echo "" >&2
    exit 1
else
input=$1
fi

# In case the user does not provide any min stag we provide a default value of 1
################################################################################
if [ ! -n "$2" ]
then
    minstag=1
else
    minstag=$2
fi

# 1) In each experiment, select the junctions supported by at least n staggered or just 1 staggered + PE information 
####################################################################################################################
echo I am selecting the junctions supported by at least n staggered or just 1 staggered + PE information, reads \in each experiment >&2
cat $input | while read lid f
do
awk -v minstag=$minstag '(($(NF-2)!=".")||($2>=minstag)){$(NF-2)=""; split($1,a,":"); split(a[1],a1,"_"); split(a[2],a2,"_"); if((a1[3]!=".")&&(a2[3]!=".")){if(($(NF-1)==".")||(($(NF-1)<80)||($NF<30))){$(NF-1)=""; $NF=""; print}}}' $f > ${f%.txt}\_$minstag\staggered_withmaxbegandend.txt
done 

# 2) Make the junctions obtained in all experiments with their most extreme beg and end across all experiments
##############################################################################################################
echo I am determining the set of junctions found \in all experiments with their most extreme beg and end \in all exp >&2
echo as well as with all other information from the input file that is not expression >&2
cat $input | while read lid f
do
cat ${f%.txt}\_$minstag\staggered_withmaxbegandend.txt
done | awk -v fld1=4 -v fld2=5 -f $BEGEND > allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_allinfo.txt

# 3) Add the number of staggered reads for each experiment together with a header to the whole matrix
#####################################################################################################
echo I am adding the number of staggered reads for each experiment >&2
awk 'BEGIN{print "junc_id beg end samechrstr okgxorder dist ss1 ss2 gnlist1 gnlist2 gnname1 gnname2 gnbt1 gnbt2"} {print $0}' allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_allinfo.txt > allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_allinfo_withexpr_eachexp.txt
cat $input | while read lid f
do
awk -v lid=$lid -v fileRef=$f 'BEGIN{while (getline < fileRef >0){nb[$1]=$2}} NR==1{print $0, lid}NR>=2{print $0, (nb[$1]!="" ? nb[$1] : 0)}' allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_allinfo_withexpr_eachexp.txt > tmp
mv tmp allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_allinfo_withexpr_eachexp.txt
done

# 4) Clean
##########
echo I am cleaning >&2
cat $input | while read lid f
do
rm ${f%.txt}\_$minstag\staggered_withmaxbegandend.txt
done
rm allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_allinfo.txt
echo I am done >&2
