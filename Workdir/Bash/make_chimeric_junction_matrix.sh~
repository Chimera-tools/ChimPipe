#!/bin/bash

# ~brodriguez/Chimeras_project/Chimeras_detection_pipeline/Chimera_mapping/Workdir/Bash/make_chimeric_junction_matrix.sh
########################################
# this script takes as input the chimeric junctions obtained from chimsplice 
# (ie from find_exon_exon_connections_from_splitmappings_better2.sh followed by 
# find_chimeric_junctions_from_exon_to_exon_connections_better2.sh) and produces
# a matrix with all the chimeric junctions obtained for all exp as well as their 
# beg and end across all experiments, whether the two parts are on the same chr 
# and strand, whether they are on the expected genomic order, their distance, 
# the list of genes which exons overlap each part of the junction, the real names 
# of those genes, the biotypes of those genes and the number of staggered reads
# supporting the junction in each exp.


# !!! should not be run twice in the same dir at the same time(common files) !!!

# More precisely, this script takes as input:
#############################################
# 1) a 3 column file containing the chimeric junction files associated to each exp
#   = where column 1 is a unique identifier for the exp, column 2 is the absolute
#   path to the chimeric junction file with the number of staggered reads, and
#   column 3 is the file of block pairs for the junctions = with the exact beg and end
#   example
#   LID20728 /users/rg/sdjebali/ENCODE_AWG/Analyses/Mouse_Human/Chimeras/Mouse/LID20728/FromTotal/distinct_junctions_and_nbstaggeredsplimappings_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt /users/rg/sdjebali/ENCODE_AWG/Analyses/Mouse_Human/Chimeras/Mouse/LID20728/FromTotal/distinct_staggered_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt
#   60 (3 fields) 
# 2) a minimum number of staggered reads for a junction to be reported in the matrix
# 3) an gff annotation file with at least the exons and the gene in column 10
# 4) 1 for stranded or 0 for unstranded data
# and this script will produce as output a matrix like this one where the script has been run
#############################################################################################
# junc_id beg end samechrstr okgxorder dist gnlist1 gnlist2 gnname1 gnname2 gnbt1 gnbt2 LID8963 LID8964 LID8965 LID8966 LID8969 LID8970 LID44498 LID44499 LID16629 LID16630 LID8461 LID8462 LID16633 LID16634 LID16635 LID16636 LID8710 LID8711 LID8463 LID8464 LID45016 LID45017 LID16627 LID16628 LID8686 LID8687 LID44594 LID44497 LID16631 LID16632 LID8692 LID8701 LID46598 LID46599 LID8967 LID8968
# chr17_5405057_-:chr17_5403432_- 5405127 5403366 1 1 1625 ENSG00000091592.8, ENSG00000170233.3, NLRP1, NLRP1, protein_coding, ambiguous_orf, 0 0 0 0 0 1 14 11 0 2 0 0 0 0 0 0 0 0 2 1 2 0 4 2 0 0 47 9 2 3 0 0 3 3 1 0


# usage
#######
# make_chimeric_junction_matrix.sh file_with_junc_each_exp.txt minstag annot stranded

# example
#########
# cd ~/ENCODE_AWG/Analyses/Mouse_Human/Chimeras/Human/AllCases
# time make_chimeric_junction_matrix.sh lid_junctionfile_blockfile.txt 2 /users/rg/projects/encode/scaling_up/whole_genome/Gencode/version10/Long/Element/gen10.long.exon.gtf 1

# Notes
#######
# - Made for using on a 64 bit linux architecture
# - uses awk scripts




# Programs
###########
CUTGFF=~brodriguez/Chimeras_project/Chimeras_detection_pipeline/Chimera_mapping/Workdir/Awk/cutgff.awk
GFF2GFF=~brodriguez/Chimeras_project/Chimeras_detection_pipeline/Chimera_mapping/Workdir/Awk/gff2gff.awk
OVERLAP=~brodriguez/Chimeras_project/Chimeras_detection_pipeline/Chimera_mapping/Workdir/bin/overlap

# In case the user does not provide any input file
###################################################
if [ ! -n "$1" ]
then
    echo "" >&2
    echo Usage: make_chimeric_junction_matrix.sh file_with_junc_each_exp.txt minstag annot stranded >&2
    echo "" >&2
    echo Example: make_chimeric_junction_matrix.sh lid_junctionfile_blockfile.txt 2 /users/rg/projects/encode/scaling_up/whole_genome/Gencode/version10/Long/Element/gen10.long.exon.gtf >&2
    echo "" >&2
    echo Takes as input:
    echo 1\) a 3 column file containing the chimeric junction files associated to each exp >&2
    echo where column 1 is a unique identifier for the exp, column 2 is the absolute >&2
    echo path to the chimeric junction file with the number of staggered reads, and >&2
    echo column 3 is the file of block pairs for the junctions = with the exact beg and end >&2
    echo 2\) a minimum number of staggered reads for a junction to be reported in the matrix >&2
    echo 3\) an gff annotation file with at least the exons and the gene in column 10 file listing >&2
    echo and produces a matrix containing the chimeric junctions obtained across all experiments >&2
    echo with a lot of additional information about their position, distance, genes and expression >&2
    echo "4\) 1 for stranded or 0 for unstranded data" >&2
    echo "" >&2
    exit 1
else
input=$1
fi



# In case the user does not provide any annotation file
########################################################
# or output dir or strandedness we provide default values
#########################################################
if [ ! -n "$2" ]
then
    minstag=1
    annot=/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version15/gencode.v15.annotation.gtf
else
    minstag=$2
    if [ ! -n "$3" ]
    then
		annot=/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version15/gencode.v15.annotation.gtf
    else
		annot=$3
		if [ ! -n "$4" ]
		then
			stranded=1
		else
			stranded=$4
    fi
fi

# 0) Make an exon file with 20 fields formatted ok for overlap and determine the fields of 
##########################################################################################
# the number of the column where gene_name and gene_type information are located
################################################################################
echo I am making the exon gff file and determining where the gene_name and gene_type lie \in this gff >&2
annotbase=`basename ${annot%.gtf}`
awk '$3=="exon"' $annot | awk -v to=20 -f $CUTGFF > $annotbase.exons.20flds.gff
str=`awk 'NR==10{fldgnname="."; fldgnbt="."; k=9; while(k<=(NF-1)){if($k=="gene_name"){fldgnname=(k+1);}else{if($k=="gene_type"){fldgnbt=(k+1);}} k+=2}}END{print fldgnname, fldgnbt}' $annotbase.exons.20flds.gff`
fldgnname=`echo $str | awk '{print $1}'` 
fldgnbt=`echo $str | awk '{print $2}'`
echo they are $fldgnname and $fldgnbt respectively >&2

# 1) In each experiment, select the junctions supported by the min number of staggered reads, 
#############################################################################################
# make a nr list of them and add their min beg and their max end in the exp
###########################################################################
echo I am selecting the junctions supported by at least $minstag staggered reads \in each experiment and finding their min beg and max end >&2
cat $input | while read lid f1 f2
do
awk -v minstag=$minstag '$2>=minstag{print $0}' $f1 | awk -v fileRef=$f2 'BEGIN{while(getline < fileRef >0){split($1,a,":"); split(a[1],a1,"_"); split(a[2],a2,"_"); if (stranded==0){c=a1[1]"_"a1[3]":"a2[1]"_"a2[2]}else{c=a1[1]"_"a1[3]"_"a1[4]":"a2[1]"_"a2[2]"_"a2[4]}; if((beg[c]=="")||(a1[2]<beg[c])){beg[c]=a1[2]} if((end[c]=="")||(a2[3]>end[c])){end[c]=a2[3]} }} {print $0, beg[$1], end[$1]}' > ${f1%.txt}\_$minstag\staggered_withmaxbegandend.txt
done

# 2) Make the junctions obtained in all experiments with their min beg and max end across exp
############################################################################################# 
echo I am determining the set of junctions found \in all experiments >&2
cat $input | while read lid f1 f2
do
cat ${f1%.txt}\_$minstag\staggered_withmaxbegandend.txt
done | awk '{if((beg[$1]=="")||($3<beg[$1])){beg[$1]=$3} if((end[$1]=="")||($4>end[$1])){end[$1]=$4}} END{for(j in beg){print j, beg[j], end[j]}}' > allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend.txt

# 3) Add 2 booleans, one for intrachrsamestr and one for ok genomic order when same intrachrsamestr 
###################################################################################################
# as well as the distance when it makes sense
#############################################
# !!! note that whatever the strand, the two parts of the junctions should always be in 5' to 3' order wrt + strand !!!
# !!! this is due to the bam convention !!!
# !!! note that this should be changed at some point !!!
echo I am adding the information of intrachrstr, okgxorder and distance when it makes sense >&2
awk '{split($1,a,":"); split(a[1],a1,"_"); split(a[2],a2,"_"); samechrstr=(a1[1]==a2[1]&&(a1[3]==a2[3]) ? 1 : 0); okgxorder=((samechrstr==1) ? (((a2[2]-a1[2])>=0) ? 1 : 0) : "NA"); print $0, samechrstr, okgxorder, ((okgxorder==1) ? (a2[2]-a1[2]) : "NA");}' allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend.txt > allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_samechrstr_okgxorder_dist.txt

# 4) add the list of gene ids, real gene names, gene biotypes for each part of each junction
############################################################################################
echo I am adding the list of gene ids, gene names and gene biotypes for each part of the junction >&2
# a. first make the gff file of the chimeric junctions
######################################################
echo I am making the gff file of the junctions >&2
awk '{split($1,a,":"); split(a[1],a1,"_"); split(a[2],a2,"_"); if (stranded==0){print a1[1], ".", ".", $2, a1[2], ".", ".", ".", "junc:", $1; print a2[1], ".", ".", a2[2], $3, ".", ".", ".", "junc:", $1}else{print a1[1], ".", ".", $2, a1[2], ".", a1[3], ".", "junc:", $1; print a2[1], ".", ".", a2[2], $3, ".", a2[3], ".", "junc:", $1}}' allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_samechrstr_okgxorder_dist.txt | awk -f $GFF2GFF > allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_two_parts.gff
# b. then overlap to know which genes have their exons overlapping each part of the junction 
#############################################################################################
echo I am overlapping the two parts of each junction with the annotated exons>&2
$OVERLAP allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_two_parts.gff $annotbase.exons.20flds.gff -st stranded -m 10 -f ex -nr -o allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_two_parts_withgnlist.gff
# c. then add the gene real names and the gene biotypes 
#######################################################
echo I am adding the gene real names and biotypes >&2
awk -v gnname=$fldgnname -v gnbt=$fldgnbt -v fileRef=$annotbase.exons.20flds.gff 'BEGIN{while (getline < fileRef >0){split($gnbt,a,"\""); split($gnname,b,"\""); bt[$10]=a[2]; name[$10]=b[2]}} {s1=""; s2=""; split($NF,a,","); k=1; while(a[k]!=""){s1=(s1)(bt[a[k]])(","); s2=(s2)(name[a[k]])(","); k++;} print $0, "gnnamelist:", s2, "btlist:", s1}' allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_two_parts_withgnlist.gff | awk -f $GFF2GFF > allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_two_parts_withgnlist_gnnamelist_btlist.gff 
# d. finally introduce all this information in the matrix
#########################################################
echo I am putting this information back to the matrix >&2
awk -v fileRef=allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_two_parts_withgnlist_gnnamelist_btlist.gff 'BEGIN{while (getline < fileRef >0){split($10,a,":"); split(a[1],a1,"_"); split(a[2],a2,"_"); s=""; split($14,b,","); k=1; while(b[k]!=""){split(b[k],c,"\""); s=(s)(c[2])(","); k++;} if($5==a1[2]){gnlist[$10,1]=s; gnnamelist[$10,1]=$16; btlist[$10,1]=$18;}else{if($4==a2[2]){gnlist[$10,2]=s; gnnamelist[$10,2]=$16; btlist[$10,2]=$18}}}} {print $0, gnlist[$1,1], gnlist[$1,2], gnnamelist[$1,1], gnnamelist[$1,2], btlist[$1,1], btlist[$1,2]}'  allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_samechrstr_okgxorder_dist.txt > allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_samechrstr_okgxorder_dist_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2.txt 

# 5) Add the number of staggered reads for each experiment together with a header to the whole matrix
#####################################################################################################
echo I am adding the number of staggered reads for each experiment >&2
awk 'BEGIN{print "junc_id beg end samechrstr okgxorder dist gnlist1 gnlist2 gnname1 gnname2 gnbt1 gnbt2"} {print $0}' allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_samechrstr_okgxorder_dist_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2.txt > allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_samechrstr_okgxorder_dist_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_withexpr_eachexp.txt
cat $input | while read lid f1 f2
do
awk -v lid=$lid -v fileRef=$f1 'BEGIN{while (getline < fileRef >0){nb[$1]=$2}} NR==1{print $0, lid}NR>=2{print $0, (nb[$1]!="" ? nb[$1] : 0)}' allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_samechrstr_okgxorder_dist_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_withexpr_eachexp.txt > tmp
mv tmp allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_samechrstr_okgxorder_dist_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_withexpr_eachexp.txt
done

# 6) Clean
##########
echo I am cleaning >&2
# rm $annotbase.exons.20flds.gff 
# rm allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend.txt 
# rm allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_two_parts.gff allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_two_parts_withgnlist.gff
# rm allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_samechrstr_okgxorder_dist.txt allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_samechrstr_okgxorder_dist_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2.txt
echo I am done >&2
