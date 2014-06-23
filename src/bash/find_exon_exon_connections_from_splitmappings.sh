#!/bin/bash

# !!! assumes /1 and /2 for pe stranded data 
# note that this script heavily depends on the format of the split mapping gff files
# and expects a format like this
# chr4	C434.1.unmapped.gem.split-map.gz	splitread	103592073	103592124	.	-	.	ID "HWI-ST858_57:2:1101:10001:72658#15@0/1"; Identity "1:0:0"; Split "[52]"
# chr4	C434.1.unmapped.gem.split-map.gz	splitread	103590183	103590206	.	-	.	ID "HWI-ST858_57:2:1101:10001:72658#15@0/1"; Identity "1:0:0"; Split "[52]"
# with 14 fields and the 2 parts of the split mapped reads following each other
# note that this script only asks the two parts of the split mapped read to (strandedly when data is stranded)
# overlap an exon, not more.
# Note: the annotation file should only have exons and the gene id must be in field no 10
# !!! be careful: in fact I realize that I could sort the file1 of overlap beforehand
# since the order will be retrieved afterwards, and then no need to specify and input 
# in overlap a sorting directory ....
# usage
#######
# find_exon_exon_connections_from_splitmappings_better2.sh split_mapping_files_paths.txt annot.gtf outputdir strandedness 

# example
#########
# cd ~/ENCODE_AWG/Analyses/ChimericConnections
# time find_exon_exon_connections_from_splitmappings_better2.sh 001N/split_mapping_files_exp_001N.txt /users/rg/projects/encode/scaling_up/whole_genome/Gencode/version15/gencode.v15.annotation.gtf 001N 1 2> 001N/find_exon_exon_connections_from_splitmappings.err

# - Made for using on a 64 bit linux architecture
# - uses awk scripts

function usage
{
cat <<instructions
	Usage:    find_exon_exon_connections_from_splitmappings.sh split_mapping_files_paths.txt annot.gtf outputdir strandedness 
    
    Example:  find_exon_exon_connections_from_splitmappings.sh 001N/split_mapping_files_exp_001N.txt /users/rg/projects/encode/scaling_up/whole_genome/Gencode/version15/gencode.v15.annotation.gtf 001N 1 
    
    Takes a file listing the absolute paths to split-mapping files \(normally .gtf.gz\) of an experiment, an annotation, 
    an output directory and the strandedness of the dataset (0 if unstranded, any positive value otherwise\) and produces
    intermediate and final files there including an exon to exon connection file for each input file 
	exit 0
instructions
}


# GETTING INPUT ARGUMENTS 
#########################
input=$1
annot=$2
outdir=$3
stranded=$4


# SETTING VARIABLES AND INPUT FILES
###################################
if [[ ! -e $input ]]; then printf "\n\tERROR: Please specify a valid input file\n\n" >&2; usage; fi
if [[ ! -e $annot ]]; then printf "\n\tERROR:Please specify a valid annotation file\n\n" >&2; usage; fi
if [[ ! -d $outdir ]]; then outdir=.; fi
if [[ $stranded == "" ]]; then stranded=0; fi

# DIRECTORIES 
#############
# IMPORTANT! rootDir is an environmental variable defined and exported in the main script "ChimPipe.sh" which contains the path to the root folder of ChimPipe pipeline. 
binDir=$rootDir/bin
awkDir=$rootDir/src/awk
bashDir=$rootDir/src/bash

# PROGRAMS
##########
CUTGFF=$awkDir/cutgff.awk
GFF2GFF=$awkDir/gff2gff.awk
OVERLAP=$binDir/overlap-3.1/overlap

# START
########

# Produces an exon file with 12 fields only from the complete annotation file
##############################################################################
echo I am making an exon file with 12 fields >&2
annotbase=`basename ${annot%.gtf}`
awk '$3=="exon"' $annot | awk -v to=12 -f $CUTGFF > $annotbase.exons.12flds.gff

# For each input file it does 7 actions (time taken is specified in comments). Input files are supposed to be .gtf.gz
#####################################################################################################################
cat $input | while read f
do
b=`basename ${f%.gz}`
btmp=${b%.gtf}

echo == For input file $b == >&2

# 0) determine the number of fields in the file
################################################
echo I am determining the number of fields of the input file >&2
nbfld=`zcat $f | head -1 | awk '{print NF}'`
echo it is $nbfld >&2
 
# 1) all split mappings = gff file that goes two lines by two lines -> 2 sec
#############################################################################
echo I am making a split-mapping file that goes 2 lines by 2 lines >&2
zcat $f -c > $outdir/$b
  
# 2.1) splitmapping part 1 = gff file in same order as split mapping file above
################################################################################
# 2.2) splitmapping part 2 = gff file in same order as split mapping file above -> 3 sec
#########################################################################################
echo I am making two gff files out of it: one for splitmapping part 1 and one for splitmapping part 2 >&2
awk -v outdir=$outdir -v btmp=$btmp 'NR%2==1{print > outdir"/"btmp"_part1.gff"}NR%2==0{print > outdir"/"btmp"_part2.gff"}' $outdir/$b


# 3.1) 2.1 with list of strandedly overlapping exons in same order as 2.1
#########################################################################
# 3.2) 2.2 with list of strandedly overlapping exons in same order as 2.2 -> 45 sec for grape unique split mappings 
###################################################################################################################
#      (36 minutes on all cshl 2 block alignments and produces huge files)
##########################################################################
echo I am adding to these files list of strandedly overlapping exons and keep the same order as \in input file >&2
$OVERLAP $outdir/$btmp\_part1.gff $annotbase.exons.12flds.gff -m -1 -st $stranded -f ex -nr -v | awk -v fileRef=$annotbase.exons.12flds.gff 'BEGIN{while (getline < fileRef >0){split($10,a,"\""); gnlist[$1"_"$4"_"$5"_"$7]=(gnlist[$1"_"$4"_"$5"_"$7])(a[2])(",");}} {s=""; if($NF=="."){s=".";} else{split($NF,a,","); k=1; while(a[k]!=""){s=(s)(gnlist[a[k]])(","); k++;}} gsub(/,,/,",",s); print $1"_"$4"_"$5"_"$7, $NF, s;}' > $outdir/$btmp\_part1_coord_exlist_gnlist.txt
$OVERLAP $outdir/$btmp\_part2.gff $annotbase.exons.12flds.gff -m -1 -st $stranded -f ex -nr -v | awk -v fileRef=$annotbase.exons.12flds.gff 'BEGIN{while (getline < fileRef >0){split($10,a,"\""); gnlist[$1"_"$4"_"$5"_"$7]=(gnlist[$1"_"$4"_"$5"_"$7])(a[2])(",");}} {s=""; if($NF=="."){s=".";} else{split($NF,a,","); k=1; while(a[k]!=""){s=(s)(gnlist[a[k]])(","); k++;}} gsub(/,,/,",",s); print $1"_"$4"_"$5"_"$7, $NF, s;}' > $outdir/$btmp\_part2_coord_exlist_gnlist.txt
awk -v fileRef=$outdir/$btmp\_part1_coord_exlist_gnlist.txt 'BEGIN{while (getline < fileRef >0){exonlist[$1]=$2; genelist[$1]=$3}} {print $0, "exlist", exonlist[$1"_"$4"_"$5"_"$7], "gnlist", genelist[$1"_"$4"_"$5"_"$7]}' $outdir/$btmp\_part1.gff > $outdir/$btmp\_part1_withexlist_gnlist.gff
awk -v fileRef=$outdir/$btmp\_part2_coord_exlist_gnlist.txt 'BEGIN{while (getline < fileRef >0){exonlist[$1]=$2; genelist[$1]=$3}} {print $0, "exlist", exonlist[$1"_"$4"_"$5"_"$7], "gnlist", genelist[$1"_"$4"_"$5"_"$7]}' $outdir/$btmp\_part2.gff > $outdir/$btmp\_part2_withexlist_gnlist.gff

# 4) 1) when both part1 and part2 overlap an exon, with indication of
#####################################################################
#    - part1 coord
#    - part2 coord
#    - readid
#    - list of exons seen by part1
#    - list of exons seen by part2  -> 2 sec for grape unique split mappings (5,5 minutes for all cshl 2 block alignments)
# chr1	hts	alblock	14755	14829	.	-	.	name: "HWI-ST935:126:C16R7ACXX:3:1311:3519:22382/2"; exlist chr1_14363_14829_-, gnlist ENSG00000227232.3,ENSG00000227232.3,ENSG00000227232.3,	chr1	hts	alblock	14970	14970	.	-	.	name: "HWI-ST935:126:C16R7ACXX:3:1311:3519:22382/2"; exlist chr1_14970_15038_-, gnlist ENSG00000227232.3,ENSG00000227232.3,ENSG00000227232.3,
echo Whenever part1 and part2 overlap an exon I am reporting their coord together with the read id and the list of such exons for each part >&2
echo Note that I am only reporting split mappings where part1 and part2 overlap exons belonging to different genes >&2
paste $outdir/$btmp\_part1_withexlist_gnlist.gff $outdir/$btmp\_part2_withexlist_gnlist.gff | awk -v nbfld=$nbfld '($(nbfld+2)!=".")&&($((nbfld+2)*2+2)!="."){ok=1; split($(nbfld+4),a,","); split($((nbfld+4)*2),b,","); k=1; while((ok==1)&&(a[k]!="")){l=1; while((ok==1)&&(b[l]!="")){if(a[k]==b[l]){ok=0} l++} k++} if(ok==1){split($10,a,"\""); print $1"_"$4"_"$5"_"$7, $(nbfld+5)"_"$(nbfld+8)"_"$(nbfld+9)"_"$(nbfld+11), a[2], $(nbfld+2), $((nbfld+2)*2+2);}}' | sort | uniq > $outdir/$btmp\_with_two_parts_overex_coord1_coord2_read_listex1_listex2.txt


# 5) from 4) make all the A -> B connections between exons where there is a split mapping which first part
##########################################################################################################
#    strandedly overlaps A and second part strandedly overlaps B and for each of them list
##########################################################################################
#    - exonA coord
#    - exonB coord
#    - list of splitmapping reads detecting this connection
#    - list of corresponding split mapping part 1
#    - list of corresponding split mapping part 2  (12 sec for grape unique split mappings and 63 minutes for all cshl 2 block alignments)
echo From the previous file I am making all the exon A -\> exon B connections where there is a splitmapping >&2
echo which first part strandedly overlaps A and second part strandedly overlaps B >&2
awk '{split($4,a,","); split($5,b,","); k=1; while(a[k]!=""){l=1; while(b[l]!=""){nb[a[k]"?"b[l]]++; rl[a[k]"?"b[l]]=(rl[a[k]"?"b[l]])($3)(","); sm1[a[k]"?"b[l]]=(sm1[a[k]"?"b[l]])($1)(","); sm2[a[k]"?"b[l]]=(sm2[a[k]"?"b[l]])($2)(","); l++;} k++;}}END{for(c in nb){split(c,c1,"?"); print c1[1], c1[2], rl[c], sm1[c], sm2[c]}}' $outdir/$btmp\_with_two_parts_overex_coord1_coord2_read_listex1_listex2.txt > $outdir/exonA_exonB_with_splitmapping_part1overA_part2overB_readlist_sm1list_sm2list_$btmp.txt

# 6) 5) with a list of staggered split-mappings and a list of the total number of splitmappings supporting the exon-exon junction defined by their coordinates in both cases   
#############################################################################################################################################################################
#    eg chr10_101948896_101948915_-:chr19_46190733_46190788_-,chr10_101948897_101948915_-:chr19_46190733_46190789_-,chr10_101948898_101948915_-:chr19_46190733_46190790_-  
#    -> 8 min and LONGEST STEP for grape unique split mappings but 13 minutes for all cshl 2 block alignments
# note: it is unclear to me why I have to split twice according to ",", will have to go step by step
# to understand it again
echo I am reporting the same info but adding the list of staggered split-mappings defined by their coordinates >&2
cat $outdir/exonA_exonB_with_splitmapping_part1overA_part2overB_readlist_sm1list_sm2list_$btmp.txt | while read ea eb rl sm1 sm2
do
staggered=`echo $sm1 $sm2 | awk '{split($1,a,","); split($2,b,","); k=1; while(a[k]!=""){split(a[k],a1,","); split(b[k],b1,","); l=1; while(a1[l]!=""){print a1[l]":"b1[l]; l++;} k++;}}' | sort | uniq | awk '{s=(s)($1)(",")}END{print s}'`
total=`echo $sm1 $sm2 | awk '{split($1,a,","); split($2,b,","); k=1; while(a[k]!=""){split(a[k],a1,","); split(b[k],b1,","); l=1; while(a1[l]!=""){print a1[l]":"b1[l]; l++;} k++;}}' | awk '{s=(s)($1)(",")}END{print s}'` 
echo $ea $eb $rl $sm1 $sm2 $staggered $total
done > $outdir/exonA_exonB_with_splitmapping_part1overA_part2overB_readlist_sm1list_sm2list_staggeredlist_totalist_$btmp.txt

# 7) cleaning and zipping -> 6 seconds
######################################
echo I am cleaning intermediate files and zipping the ones that are not of immediate use >&2

rm $outdir/$b $outdir/$btmp\_part1.gff $outdir/$btmp\_part2.gff $outdir/$btmp\_part1_coord_exlist_gnlist.txt $outdir/$btmp\_part2_coord_exlist_gnlist.txt $outdir/$btmp\_part1_withexlist_gnlist.gff $outdir/$btmp\_part2_withexlist_gnlist.gff $outdir/$btmp\_with_two_parts_overex_coord1_coord2_read_listex1_listex2.txt $outdir/exonA_exonB_with_splitmapping_part1overA_part2overB_readlist_sm1list_sm2list_$btmp.txt 
gzip -f $outdir/exonA_exonB_with_splitmapping_part1overA_part2overB_readlist_sm1list_sm2list_staggeredlist_totalist_$btmp.txt 

done

# more cleaning
rm $annotbase.exons.12flds.gff 

echo I am finished >&2


