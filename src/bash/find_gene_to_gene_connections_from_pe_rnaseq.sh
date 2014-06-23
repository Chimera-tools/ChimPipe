#!/bin/bash

# ~/bin/find_gene_to_gene_connections_from_pe_rnaseq_fast.sh
# same aim as find_gene_to_gene_connections_from_pe_rnaseq_gal.sh but using bedtools
# intersectBed instead of overlap to be faster. Another important difference is that
# we use all the exons from the annotation file provided, not only the exons belonging
# to one gene only. A difference from the I/O point of view is the number and order of 
# the arguments, that is different from before and that is now the same as for the
# 2 subscripts of chimsplice, except that instead of a txt file with the gff files
# we now directly take a bam file. A final difference is that here the mappings could
# just (strandedly if stranded data) overlap the exon, and not be totally included in it
# like in find_gene_to_gene_connections_from_pe_rnaseq_gal.sh

# usage:
########
# find_gene_to_gene_connections_from_pe_rnaseq_fast.sh mapping.bam annot.gff outputdir strandedness
# where:
# - mapping.bam is the pe mapping file (compulsory)
# - annot.gff is the annotation file, which should at least have exons (default v15)
# - outputdir is the directory in which the user wants the results to be stored (default cwd)
# - strandedness is a boolean saying whether the data is stranded (default not)
# NOTE that the mapping must be paired end and the read names must end by /1 and /2 (should be more general in the future)
# NOTE that we consider all exons of the annotation and the mappings only have to overlap them by 1 bp (not inclusion)
# NOTE assumes that the gene id in the gff file is in column 10
# NOTE assumes that the geneid does not contain a dash


# In case the user does not provide any input file, an error message is raised
##############################################################################
if [ ! -n "$1" ]
then
    echo "" >&2
    echo Usage:    find_gene_to_gene_connections_from_pe_rnaseq_fast.sh mapping.bam annot.gff outputdir strandedness >&2
    echo "" >&2
    echo Example:  find_gene_to_gene_connections_from_pe_rnaseq_fast.sh /users/rg/dgonzalez/Projects/TripleMama/hg19/TM05B9085PE/SAM/05B9085.merged.bam /users/rg/dgonzalez/ReferenceAnnotations/H.sapiens/gencode.v6.annot.plus.tRNAs.female.gtf ~/BreastCancer/Selection_for_RTPCR/Confirm_by_Solexa/05B9085_map.PE.40 0 >&2
    echo "" >&2
    echo Takes a mapping file \in bam format of an experiment, an annotation \in gff format, >&2
    echo "an output directory and the strandedness of the dataset \(0 if unstranded, any positive value otherwise\)" >&2
    echo "and produces the file of \(directed\) gene to gene connections found by the pe mappings with the number" >&2
    echo "of mappings supporting the connection\." >&2 
    echo "For a connection g1 to g2 to exist there must be at least one mapping where the first mate is \(strandedly" >&2
    echo "if data is stranded\) overlapping with an exon of g1 and the second mate is \(strandedly if data is stranded\)" >&2
    echo overlapping with an exon of g2 >&2
    echo "" >&2
    exit 0
fi

bamfile=$1

# In case the user does not provide any annotation file
########################################################
# or output directory or strandedness, default values are provided
###################################################################
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

# Directories
#############
# IMPORTANT! rootDir is an environmental variable defined and exported in the main script "ChimPipe.sh" which contains the path to the root folder of ChimPipe pipeline. 
awkDir=$rootDir/src/awk
binDir=$rootDir/bin

# Programs
##########
BEDTOOLS=$binDir/bedtools2-2.20.1/bin/bedtools
CUTGFF=$awkDir/cutgff.awk
REMOVEREDUND=$awkDir/remove_redund.awk
GFF2GFF=$awkDir/gff2gff.awk

# Intersect the bam file with the exons of the annotation (longest step)
########################################################################
# treat the split mappings as separate blocks
#############################################
echo I am intersecting the bam file with the exons of the annotation and providing >&2
echo for each mapping of a read the read id with the non redundant list of genes whose >&2
echo exons are overlapped by the read \(longest step, overlap done strandedly if data is stranded\) >&2
if [ $stranded == 0 ]
then
str=""
else
str="-s"
fi
awk '$3=="exon"' $annot | awk -v to=12 -f $CUTGFF | $BEDTOOLS intersect -abam $bamfile -b stdin $str -split -bed -wao | awk '$NF>0{split($22,a,"\""); gnlist[$4]=(gnlist[$4])(a[2])(",");}END{for(m in gnlist){print m, gnlist[m]}}' | awk -v fldlist=gnlist:2 -f $REMOVEREDUND | awk '{split($4,a,","); k=1; s=""; while(a[k]!=""){split(a[k],b,":"); s=(s)(b[1])(","); k++} print $1, s}' | gzip > $outdir/readid_gnlist_whoseexoverread_noredund.txt.gz
echo done >&2

# For each read that has both mates with a gene list, provide all the pairs 
###########################################################################
# of different genes where the first one is in the gene list of /1 and the 
#########################################################################
# second one is in the gene list of /2
#######################################
echo For each read that has both mates associated to a gene list >&2
echo I am making all the pairs of different genes where the first one is>&2
echo \in the list of the first mate and where the second one >&2
echo is \in the list of the second mate >&2
zcat $outdir/readid_gnlist_whoseexoverread_noredund.txt.gz | awk '{split($1,a,"/"); exists[a[1]]++; gnlist[a[1],a[2]]=$2}END{for(r in exists){s=""; if((gnlist[r,1]!="")&&(gnlist[r,2]!="")){split(gnlist[r,1],a,","); split(gnlist[r,2],b,","); k=1; while(a[k]!=""){l=1; while(b[l]!=""){if(a[k]!=b[l]){s=(s)(a[k]"-"b[l])(",");} l++} k++}} print r, s}}' | gzip > $outdir/readid_twomateswithgnlist_alldiffgnpairs_where_1stassociatedto1stmate_and2ndto2ndmate.txt.gz
echo done >&2

# Gather this information in order to report for each pair of different genes (in alphabetical order)
####################################################################################################
# the number of reads where the two mates support the pair
##########################################################
echo I am gathering this information \in order to report all the pairs of diffent genes >&2
echo \in alphabetical order\, that are supported by pe reads together with the number of reads >&2
echo supporting them >&2
zcat $outdir/readid_twomateswithgnlist_alldiffgnpairs_where_1stassociatedto1stmate_and2ndto2ndmate.txt.gz | awk '{split($2,a,","); k=1; while(a[k]!=""){split(a[k],b,"-"); if(b[1]<b[2]){print b[1], b[2]}else{print b[2], b[1]} k++}}' | sort | uniq -c | awk '{print $2, $3, $1}' > $outdir/pairs_of_diff_gn_supported_by_pereads_nbpereads.txt
echo done >&2


