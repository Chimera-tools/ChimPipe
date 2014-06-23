#!/bin/bash

# usage
#######
# find_chimeric_junctions_from_exon_to_exon_connections_better2.sh split_mapping_files_paths.txt annot.gtf outputdir strandedness

# to be used following the use of find_exon_exon_connections_from_splitmappings_better2.sh with same arguments 

# !!! be careful: cannot be run at the same time in the same directory since writes annot file with same name

# example
# cd ~/ENCODE_AWG/Analyses/ChimericConnections
# time find_chimeric_junctions_from_exon_to_exon_connections_better2.sh 001N/split_mapping_files_exp_001N.txt /users/rg/projects/encode/scaling_up/whole_genome/Gencode/version7/gencode.v7.annotation_exons.gtf 001N 1 > 001N/chimeric_junctions_report.txt 2> 001N/find_chimeric_junctions_from_exon_to_exon_connections.err
# real	0m25.021s

# produces a report like this
#############################
# total exonA to exonB connections: 323862
# number of exon A to exon B connections where A and B are in one gene respectively and those genes are different: 114332
# total number of distinct chimeric junctions: 70778
# number of distinct chimeric junctions seen by at least 10 staggered split-mappings: 243

# Directories
#############
# IMPORTANT! rootDir is an environmental variable defined and exported in the main script "ChimPipe.sh" which contains the path to the root folder of ChimPipe pipeline. 
awkDir=$rootDir/src/awk
bashDir=$rootDir/src/bash

# Programs
###########
GFF2GFF=$awkDir/gff2gff.awk
STAG=$awkDir/staggered_to_junction.awk
MATRIX=$bashDir/make_chimeric_junction_matrix.sh
SPLICE_SITES=$bashDir/find_consensus_splice_sites.sh

# Temporary directory for the sorting
######################################
tmpdir=~brodriguez/Tmp

# In case the user does not provide any input file
###################################################
if [ ! -n "$1" ]
then
    echo "" >&2
    echo Usage:    find_chimeric_junctions_from_exon_to_exon_connections.sh split_mapping_files_paths.txt annot.gtf outputdir [strandedness] >&2
    echo "" >&2
    echo Example:  find_chimeric_junctions_from_exon_to_exon_connections.sh 001N/split_mapping_files_exp_001N.txt  /users/rg/projects/encode/scaling_up/whole_genome/Gencode/version7/gencode.v7.annotation_exons.gtf 001N 1 >&2
    echo "" >&2
    echo Takes a file listing the absolute paths to split-mapping files \(normally .gtf.gz\) of an experiment, an output >&2
    echo directory and the strandedness of the data, and produces a summary statistics file as output with the number >&2
    echo of chimeric junctions seen >&2
    echo "" >&2
    exit 0
fi

# In case the user does not provide any annotation file
########################################################
# or output dir or strandedness we provide default values
#########################################################
if [ ! -n "$2" ]
then
annot=/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version15/gencode.v15.annotation.gtf
outdir=.
stranded=0
else
annot=$2
if [ ! -n "$3" ]
then
outdir=.
stranded=0
else
outdir=$3
if [ ! -n "$4" ]
then
stranded=0
else
stranded=$4
fi
fi
fi

# Exons that are in one gene only
#################################
echo I am making the file of exons that are \in one gene only from the annotation >&2
awk '$3=="exon"{seen[$1"_"$4"_"$5"_"$7,$10]++; if(seen[$1"_"$4"_"$5"_"$7,$10]==1){split($10,a,"\""); nb[$1"_"$4"_"$5"_"$7]++; genelist[$1"_"$4"_"$5"_"$7]=(genelist[$1"_"$4"_"$5"_"$7])(a[2])(",");}} END{for(e in nb){if(nb[e]==1){split(genelist[e],a,","); print e, a[1]}}}' $annot > exoncoord_inonegene_gene.txt

# Number of exon to exon connections detected by all kinds of mappings of the experiment
#########################################################################################
echo I am computing the total number of exon A to exon B connections >&2
cat $1 | while read f
do
b=`basename ${f%.gz}`
btmp=${b%.gtf}
zcat $outdir/exonA_exonB_with_splitmapping_part1overA_part2overB_readlist_sm1list_sm2list_staggeredlist_totalist_$btmp.txt.gz | awk '{print $1, $2}' | sort -T $tmpdir
done | sort -T $tmpdir | uniq | wc -l | awk '{print "total exonA to exonB connections:", $1}' 


# List the staggered splitmappings and the total splitmappings detecting each A-> B connection -> 8 sec
###########################################################################
echo I am listing the split\-mappings detecting each A\-\> B connection >&2
cat $1 | while read f
do
b=`basename ${f%.gz}`
btmp=${b%.gtf}
zcat $outdir/exonA_exonB_with_splitmapping_part1overA_part2overB_readlist_sm1list_sm2list_staggeredlist_totalist_$btmp.txt.gz
done | awk '{staggered[$1":"$2]=(staggered[$1":"$2])($6)(","); total[$1":"$2]=(total[$1":"$2])($7)(",")};  END{for(i in staggered){split(i,a,":"); gsub(/,,/,",",staggered[i]); gsub(/,,/,",",total[i]); print a[1], a[2], staggered[i], total[i]}}' > $outdir/exonA_exonB_with_splitmapping_part1overA_part2overB_staggeredlist_totalist_total.txt


# From the list of exon A -> exon B connections select the ones where both A and B are 
#######################################################################################
# in only 1 gene and where those two genes are different
########################################################
echo I am selecting the exon A to exon B connections where A and B are \in one gene respectively and those genes are different >&2
awk -v fileRef=exoncoord_inonegene_gene.txt 'BEGIN{while (getline < fileRef >0){gene[$1]=$2}} ((gene[$1]!="")&&(gene[$2]!="")&&(gene[$1]!=gene[$2])){print}' $outdir/exonA_exonB_with_splitmapping_part1overA_part2overB_staggeredlist_totalist_total.txt > $outdir/exonA_exonB_with_splitmapping_part1overA_part2overB_staggeredlist_totalist_total_only_A_B_indiffgn_and_inonegn.txt
wc -l $outdir/exonA_exonB_with_splitmapping_part1overA_part2overB_staggeredlist_totalist_total_only_A_B_indiffgn_and_inonegn.txt | awk '{print "number of exon A to exon B connections where A and B are in one gene respectively and those genes are different:", $1}'

# Compute the distinct chimeric junctions found by the staggered split-mappings, the number of split-mapping reads and the number of total reads supporting the junction,  
#########################################################################################################################################################################
# the beginning and end of the junctions. Here different according to strandedness
##################################################################################
echo I am computing the distinct chimeric junctions found by the staggered split-mappings >&2

# first the distinct staggered directed splitmappings 
#####################################################
echo first the distinct staggered directed splitmappings >&2
awk '{split($3,a,","); k=1; while(a[k]!=""){print a[k]; k++;}}' $outdir/exonA_exonB_with_splitmapping_part1overA_part2overB_staggeredlist_totalist_total_only_A_B_indiffgn_and_inonegn.txt | sort -T $tmpdir | uniq -c | awk '{print $2,$1}' > $outdir/distinct_staggered_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt

# then the distinct junctions (keep redundancy value and direction) as well as their beg and end
################################################################################################
echo then the distinct junctions \(keeping redundancy value and direction\) as their beg and end>&2
awk -v stranded=$stranded -f $STAG $outdir/distinct_staggered_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt > $outdir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt

wc -l $outdir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt | awk '{print "total number of distinct chimeric junctions:", $1}' 

# Complete the chimeric junctions matrix with additional information: whether the two parts are on the same chrand strand, whether they are on the expected genomic order, 
##########################################################################################################################################################################
# their distance, the list of genes which exons overlap each part of the junction, the real names of those genes, the biotypes of those genes and the number of staggered reads
############################################################################################################################################################################## 
# supporting the junction 
#########################
# Output: $outdir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt
echo I am completing the chimeric junctions matrix with additional information: whether the two parts are on the same chrand strand, whether they are on the expected genomic order, their distance, the list of genes which exons overlap each part of the junction, the real names of those genes, the biotypes of those genes and the number of staggered reads supporting the junction.>&2
bash $MATRIX $outdir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt  $annot $stranded $outdir

# Find canonical splice sites GT-AG or GC-AG and their complementary reverse AC-CT or GC-CT. Add this information to the matrix and for unstranded data writes the junctions
############################################################################################################################################################################
# in biological order (first part donor and second part acceptor) and add strand information to each part of the junction whenever it is possible 
####################################################################################################################################################
# Output: $outdir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt 
echo I am finding the canonical splice sites GT-AG or GC-AG, adding this information to the matrix and also for unstranded data writing the junctions in biological order and adding the strand to each part of the junction whenever it is possible >&2
bash $SPLICE_SITES $outdir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt $stranded $outdir


# Output the distinct chimeric junctions found by more than 10 staggered split-mappings
#######################################################################################
echo I am outputting the distinct chimeric junctions found by more than 10 staggered split-mappings >&2
awk '$2>=10' $outdir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt > $outdir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn_morethan10staggered.txt
wc -l $outdir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn_morethan10staggered.txt | awk '{print "number of distinct chimeric junctions seen by at least 10 staggered split-mappings:", $1}' 

#rm exoncoord_inonegene_gene.txt 
