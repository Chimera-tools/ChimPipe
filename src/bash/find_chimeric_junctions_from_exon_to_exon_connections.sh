#!/bin/bash

<<authors
*****************************************************************************
	
	find_chimeric_junctions_from_exon_to_exon_connections.sh
	
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

# will exit if there is an error or in a pipe
set -e -o pipefail

# In case the user does not provide any input file
###################################################
if [ ! -e "$1" ] 
then
    echo "" >&2
    echo Usage:    find_chimeric_junctions_from_exon_to_exon_connections.sh split_mapping_files_paths.txt genome_index.gem annot.gtf outputdir [strandedness] >&2
    echo "" >&2
    echo Example:  find_chimeric_junctions_from_exon_to_exon_connections.sh 001N/split_mapping_files_exp_001N.txt /users/rg/sdjebali/ENCODE_AWG/Analyses/Mouse_Human/Chimeras/Human/Homo_sapiens.GRCh37.chromosomes.chr.gem /users/rg/projects/encode/scaling_up/whole_genome/Gencode/version7/gencode.v7.annotation_exons.gtf 001N 1 >&2
    echo "" >&2
    echo Takes a file listing the absolute paths to split-mapping files \(normally .gtf.gz\) of an experiment, the indexed genome and >&2
    echo the reference annotation, an output directory and the strandedness of the data, and produces a summary statistics file as output >&2
    echo with the number of chimeric junctions seen >&2
    echo "" >&2
    exit -1
fi

# In case the user does not provide any indexed genome, annotation 
##################################################################
# file or output dir or strandedness we provide default values
##############################################################
if [ ! -n "$2" ]
then
	genome=/users/rg/sdjebali/ENCODE_AWG/Analyses/Mouse_Human/Chimeras/Human/Homo_sapiens.GRCh37.chromosomes.chr.gem
	annot=/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version15/gencode.v15.annotation.gtf
	outdir=.
	stranded=0
	spliceSites='(GT,AG),(GC,AG),(ATATC,A.),(GTATC,AT)'
else
	genome=$2
	if [ ! -n "$3" ]
	then
		annot=/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version15/gencode.v15.annotation.gtf
		outdir=.
		stranded=0
		spliceSites='(GT,AG),(GC,AG),(ATATC,A.),(GTATC,AT)'
	else
		annot=$3
		if [ ! -n "$4" ]
		then
			outdir=.
			stranded=0
			spliceSites='(GT,AG),(GC,AG),(ATATC,A.),(GTATC,AT)'
		else
			outdir=$4
			if [ ! -n "$5" ]
			then
				stranded=0
				spliceSites='(GT,AG),(GC,AG),(ATATC,A.),(GTATC,AT)'
			else
				stranded=$5
				if [ ! -n "$6" ]
				then
					spliceSites='(GT,AG),(GC,AG),(ATATC,A.),(GTATC,AT)'
				else
					spliceSites=$6
				fi
			fi
		fi
	fi
fi

# Directories 
#############
rootDir=/nfs/users/rg/brodriguez/Chimeras_project/Chimeras_detection_pipeline/ChimPipe
awkDir=$rootDir/src/awk
bashDir=$rootDir/src/bash

# Programs
###########
GFF2GFF=$awkDir/gff2gff.awk
MAKE_JUNCTIONS=$bashDir/make_chimjunctions.sh
MATRIX=$bashDir/make_chimeric_junction_matrix.sh

# Exons and a list of overlapping genes
#############################################
echo I am making a file of exons with the list of their overlapping genes from the annotation >&2

awk '$3=="exon"{seen[$1"_"$4"_"$5"_"$7,$10]++; if(seen[$1"_"$4"_"$5"_"$7,$10]==1){split($10,a,"\""); gnlist[$1"_"$4"_"$5"_"$7]=(gnlist[$1"_"$4"_"$5"_"$7])(a[2])(",")}}END{for(exon in gnlist){print exon, gnlist[exon]}}' $annot > $outdir/exoncoord_gnlist.txt

# Number of exon to exon connections detected by all kinds of mappings of the experiment
#########################################################################################
echo I am computing the total number of exon A to exon B connections >&2
cat $1 | while read f
do
b=`basename ${f%.gz}`
btmp=${b%.gff}

zcat $outdir/exonA_exonB_with_splitmapping_part1overA_part2overB_readlist_sm1list_sm2list_staggeredlist_totalist_$btmp.txt.gz | awk '{print $1, $2}' | sort -T $TMPDIR
done | sort -T $TMPDIR | uniq | wc -l | awk '{print "total exonA to exonB connections:", $1}' 

# List the staggered splitmappings and the total split-mappings detecting each A-> B connection -> 8 sec
########################################################################################################
echo I am listing the split\-mappings detecting each A\-\> B connection >&2
cat $1 | while read f
do
b=`basename ${f%.gz}`
btmp=${b%.gff}

zcat $outdir/exonA_exonB_with_splitmapping_part1overA_part2overB_readlist_sm1list_sm2list_staggeredlist_totalist_$btmp.txt.gz
done | awk '{staggered[$1":"$2]=(staggered[$1":"$2])($6)(","); total[$1":"$2]=(total[$1":"$2])($7)(",")};  END{for(i in staggered){split(i,a,":"); gsub(/,,/,",",staggered[i]); gsub(/,,/,",",total[i]); print a[1], a[2], staggered[i], total[i]}}' > $outdir/exonA_exonB_with_splitmapping_part1overA_part2overB_staggeredlist_totalist_total.txt

# From the list of exon A -> exon B connections select the ones where both A and B are 
#######################################################################################
# in different genes
####################
echo I am selecting the exon A to exon B connections where A and B are \in different genes >&2
awk -v fileRef=$outdir/exoncoord_gnlist.txt 'BEGIN{while (getline < fileRef >0){gnlist[$1]=$2}} {intersection=0;split(gnlist[$1],gnlistA,",");split(gnlist[$2],gnlistB,","); for (gnA in gnlistA){for (gnB in gnlistB){if((gnlistA[gnA]!="")&&(gnlistB[gnB]!="")&&(gnlistA[gnA]==gnlistB[gnB])){intersection=1}}}; if (intersection!=1){print}}' $outdir/exonA_exonB_with_splitmapping_part1overA_part2overB_staggeredlist_totalist_total.txt > $outdir/exonA_exonB_with_splitmapping_part1overA_part2overB_staggeredlist_totalist_total_only_A_B_indiffgn.txt
wc -l $outdir/exonA_exonB_with_splitmapping_part1overA_part2overB_staggeredlist_totalist_total_only_A_B_indiffgn.txt | awk '{print "number of exon A to exon B connections where A and B are in different genes:", $1}'

# Compute the distinct chimeric junctions found by the staggered split-mappings, the number of split-mapping reads and the number of total reads supporting the junction,  
#########################################################################################################################################################################
# the beginning and end of the junctions. Here different according to strandedness
##################################################################################
echo I am computing the distinct chimeric junctions found by the staggered split-mappings >&2

# first the distinct staggered directed splitmappings 
#####################################################
echo first the distinct staggered directed splitmappings >&2
awk '{split($3,a,","); k=1; while(a[k]!=""){print a[k]; k++;}}' $outdir/exonA_exonB_with_splitmapping_part1overA_part2overB_staggeredlist_totalist_total_only_A_B_indiffgn.txt | sort -T $TMPDIR | uniq -c | awk '{print $2, $1}' > $outdir/distinct_staggered_split_mappings_part1overA_part2overB_only_A_B_indiffgn.txt

# then the distinct junctions (keep redundancy value and direction) as well as their beg and end
################################################################################################
echo then the distinct junctions \(keeping redundancy value and direction\) as their consensus splice sites, maximum beginning and ending>&2
bash $MAKE_JUNCTIONS $outdir/distinct_staggered_split_mappings_part1overA_part2overB_only_A_B_indiffgn.txt $genome $stranded $spliceSites $outdir > $outdir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_ss1_ss2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt 

wc -l $outdir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_ss1_ss2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt | awk '{print "total number of distinct chimeric junctions:", $1}' 

# Complete the chimeric junctions matrix with additional information: whether the two parts are on the same chr and strand, whether they are on the expected genomic order, 
##########################################################################################################################################################################
# their distance, the list of genes which exons overlap each part of the junction, the real names of those genes, the biotypes of those genes and the number of staggered reads
############################################################################################################################################################################## 
# supporting the junction 
#########################
# Output: $outdir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt

echo I am completing the chimeric junctions matrix with additional information: whether the two parts are on the same chrand strand, whether they are on the expected genomic order, their distance, the list of genes which exons overlap each part of the junction, the real names of those genes, the biotypes of those genes and the number of staggered reads supporting the junction.>&2
bash $MATRIX $outdir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_ss1_ss2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt $annot $stranded $outdir > $outdir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt

