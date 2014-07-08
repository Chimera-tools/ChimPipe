#!/bin/bash

# will exit if there is an error or in a pipe
set -e -o pipefail

function usage
{
cat <<help
	
*** ChimPipe version $version ***

Execute ChimPipe (from paired-end RNA-Seq reads to chimeric junctions) on one RNA-Seq dataset (sample).
	
USAGE: ChimPipe.sh -i <fastq_file> -g <genome_index> -a <annotation> -q <quality> -e <sample_id> [OPTIONS]
 
IMPORTANT: By default runs in unstranded mode. If you have stranded data use the \"-s\" flag (see [OPTIONS]).

** Mandatory arguments:
	-i|--input			<INPUT_FILE>	First mate sequencing reads in FASTQ format. Pipeline designed to deal with paired-end
	 						data. Please make sure the second mate file is in the same directory as the first mate 
	 						file, and the files are named "YourSampleId_1.fastq.gz" and "YourSampleId_2.fastq.gz" 
	 						respectively. YourSampleId has to be provided with the -e argument.
	-g|--genome-index		<GEM>		Index for the reference genome in GEM format.
	-a|--annotation			<GTF>		Reference gene annotation file in GTF format.
	-q|--quality			<NUMBER>	Quality offset of the FASTQ files [33 | 64 | ignore].
	-e|--sample-id			<STRING>	Sample identifier (the output files will be named according to this id).    
	
** [OPTIONS] can be:
Reads information:
	-s|--stranded			<FLAG>		Flag to specify whether data is stranded. Default false (unstranded).
	--read-directionality		<STRING>	Directionality of the reads [MATE1_SENSE | MATE2_SENSE | MATE_STRAND_CSHL | SENSE | ANTISENSE | NONE]. Default NONE.
	--max-read-length		<NUMBER>	Maximum read length. This is used to create the de-novo transcriptome and acts as an upper bound. Default 150.
	
Mapping parameters:
	-M|--mism-contiguous-map	<NUMBER>	Maximum number of mismatches for the contiguous mapping steps with the GEM mapper. Default 4?.
	-m|mism-split-map		<NUMBER>	Maximum number of mismatches for the segmental mapping steps with the GEM rna-mapper. Default 4?.	
	--min-split-size		<NUMBER>	Minimum split size for the segmental mapping steps. Default 15.
	--stats				<FLAG>		Enable mapping statistics. Default disabled.
	
Chimeric junctions filter:
	--filter-chimeras		<STRING>	Configuration for the filtering module. Quoted string with 4 numbers separated by commas and ended in semicolom, 
							i.e. "1,2,75,50;", where:
											
								1st: minimum number of staggered reads spanning the chimeric junction.
								2nd: maximum number of paired-end reads encompassing the chimeric junction.		
								3rd: maximum similarity between the connected genes.
								4rd: minimum length of the high similar region between the connected genes.
	
							All these conditions have to be fulfilled for a chimeric junction to pass the filter. It is also possible to make 
							complex condifions by setting two different conditions where at least one of them has to be fulfilled. 
							I.e "10,0,00,00;1,1,00,00;". Default "5,0,80,30;1,1,80,30;".	
	--similarity-gene-pairs	<TEXT>			Text file containing similarity information between the gene pairs in the annotation. Needed for the filtering module 
							to discard junctions connecting highly similar genes. If not provided the junctions will not be filtered according this criteria. 

General:
	-o|--output-dir			<PATH>		Output directory. Default current working directory.
	--tmp-dir			<PATH>		Temporary directory. Default /tmp.
	-t|--threads			<PATH>		Number of threads to use. Default 1.
	-l|--log			<PATH>		Log level [error |warn | info | debug). Default info.
	--dry				<FLAG>		Test the pipeline. Writes the command to the standard output.
	--help				<FLAG>		Display usage information.

help
}

function printHeader {
    string=$1
    echo "`date` ***** $string *****"
}

function log {
    string=$1
    label=$2
    if [[ ! $ECHO ]];then
        if [[ "$label" != "" ]];then
            printf "[$label] $string"
        else
            printf "$string"
        fi
    fi
}

function run {
    command=($1)
    if [[ $2 ]];then
         ${2}${command[@]}
    else
        eval ${command[@]}
    fi
}

# Function 5. 
#############
function getoptions {
ARGS=`$getopt -o "i:g:a:q:e:sM:m:o:t:l:h" -l "input:,genome-index:,annotation:,quality:,sample-id:,stranded,mismatches-contiguous-mapping:,mismatches-split-mapping:,output-dir:,threads:,log:,help,read-directionality:,max-read-length:,max-read-length:,consensus-splice-sites:,min-split-size:,filter-chimeras:,similarity-gene-pairs:,stats,tmp-dir:,dry" \
      -n "$0" -- "$@"`

#Bad arguments
if [ $? -ne 0 ];
then
  exit 1
fi

# A little magic
eval set -- "$ARGS"

while true;
do
  case "$1" in
    -i|--input)
      if [ -n "$2" ];
      then
        input=$2
      fi
      shift 2;;

    -g|--genome-index)
      if [ -n "$2" ];
      then
        index=$2
      fi
      shift 2;;

    -a|--annotation)
      if [ -n "$2" ];
      then
        annot=$2
      fi
      shift 2;;
  
     -q|--quality)
      if [ -n $2 ];
      then
        quality=$2
      fi
      shift 2;;
    
    -e|--sample-id)
       if [ -n "$2" ];
       then
         lid=$2
       fi
       shift 2;;
    
    -s|--stranded)
       stranded=1
       shift;;
     
    --read-directionality)
      if [ -n $2 ];
      then
        readDirectionality=$2
      fi
      shift 2;;
     
    --max-read-length)
      if [ -n $2 ];
      then
        maxReadLength=$2
      fi
      shift 2;;

	 -M|--mismatches-contiguous-mapping)
       if [ -n "$2" ];
       then
         mism=$2
       fi
       shift 2;;
	
	-m|--mismatches-split-mapping)
       if [ -n "$2" ];
       then
         mismSplit=$2
       fi
       shift 2;;

	--consensus-splice-sites)
       if [ -n "$2" ];
       then
         spliceSites=$2
       fi
       shift 2;;
       
    --min-split-size)
       if [ -n "$2" ];
       then
         splitSize=$2
       fi
       shift 2;;   

	--stats)
      mapStats="-m"  
      shift ;;

    --filter-chimeras)
      if [ -n $2 ];
      then
       filterConf="$2"
      fi
      shift 2;;
      
    --similarity-gene-pairs)
      if [ -n $2 ];
      then
        simGnPairs="$2"
      fi
      shift 2;;

	-o|--output-dir)
      if [ -n $2 ];
      then
        outDir=$2
      fi
      shift 2;;
      
    --tmp-dir)
      if [ -n $2 ];
      then
        TMPDIR=$2
      fi
      shift 2;;
 
    -t|--threads)
      if [ -n $2 ];
      then
        threads=$2
      fi
      shift 2;;
 
    -l|--log)
      if [ -n $2 ];
      then
        logLevel=$2
      fi
      shift 2;;
 
    --dry)
      ECHO="echo "
      shift;;
    
    -h|--help)
      usage
      shift;;
    
    --)
      shift
      break;;
  esac
done
}

# SETTING UP THE ENVIRONMENT
############################

# ChimPipe version 
version=V0.7.0

# 1. ChimPipe's root directory
##############################
# It will be exported as an environmental variable since it will be used by every ChimPipe's scripts 
# to set the path to the bin, awk and bash directories. 
root=/nfs/users/rg/brodriguez/Chimeras_project/Chimeras_detection_pipeline/ChimPipe

export rootDir=$root 

# 2. Parse input arguments with getopt  
######################################
getopt=$root/bin/getopt-1.1.5/getopt

getoptions $0 $@ # call Function 5 and passing two parameters (name of the script and command used to call it)

# 3. Check input variables 
##########################
if [[ ! -e $input ]]; then log "Please specify a valid input file\n" "ERROR" >&2; usage; exit -1; fi
if [[ `basename ${input##*_}` != "1.fastq.gz" ]]; then log "Please check that the name of your FASTQ file ends with \"_1.fastq.gz\"\n" "ERROR" >&2; usage; exit -1; fi
if [[ ! -e $index ]]; then log "Please specify a valid genome index file\n" "ERROR" >&2; usage; exit -1; fi
if [[ ! -e $annot ]]; then log "Please specify a valid annotation file\n" "ERROR" >&2; usage; exit -1; fi
if [[ "$quality" == "" ]]; then log "Please specify the quality\n" "ERROR" >&2; usage; exit -1; fi
if [[ "$lid" == "" ]]; then log "Please specify the sample identifier\n" "ERROR" >&2; usage; exit -1; fi
if [[ "$bam" == "" ]]; then bam=0; fi
if [[ "$stranded" == "" ]]; then stranded=0; fi
if [[ "$readDirectionality" == "" ]]; then readDirectionality='NONE'; fi
if [[ "$mism" == "" ]]; then mism='4'; fi
if [[ "$maxReadLength" == "" ]]; then maxReadLength='150'; fi
if [[ "$spliceSites" == "" ]]; then spliceSites='GT+AG' ; fi
if [[ "$splitSize" == "" ]]; then splitSize='15'; fi
if [[ ! -d "$outDir" ]]; then outDir=${SGE_O_WORKDIR-$PWD}; fi
if [[ "$logLevel" == "" ]]; then logLevel='info'; fi
if [[ "$threads" == "" ]]; then threads='1'; fi
if [[ "$TMPDIR" == "" ]]; then TMPDIR='/tmp'; fi

if [[ "$filterConf" == "" ]]; then 			# Checking filtering module configuration
	filterConf="5,0,80,30;1,1,80,30;";		# Default
else
	if [[ ! "$filterConf" =~ ^([0-9]+,[0-9]+,[0-9]{,3},[0-9]+;){1,2}$ ]]; then
		log "Please check your filtering module configuration. Option -F\n\n" "ERROR" >&2; 
		usage; 
		exit -1;
	fi
fi

# 4. Directories
################
binDir=$rootDir/bin
awkDir=$rootDir/src/awk
bashDir=$rootDir/src/bash
if [[ ! -d $outDir/SecondMapping ]]; then mkdir $outDir/SecondMapping; fi
if [[ ! -d $outDir/FromFirstBam ]]; then mkdir $outDir/FromFirstBam; fi
if [[ ! -d $outDir/FromSecondMapping ]]; then mkdir $outDir/FromSecondMapping; fi
if [[ ! -d $outDir/Chimsplice ]]; then mkdir $outDir/Chimsplice; fi
if [[ ! -d $outDir/PE ]]; then mkdir $outDir/PE; fi

# The temporary directory will be exported as an environmental variable since it will 
# be used by every ChimPipe's scripts 
export TMPDIR=$TMPDIR

# 5. Programs/Scripts
##################
# Bash 
pipeline=$bashDir/first_mapping_pipeline.sh
chim1=$bashDir/find_exon_exon_connections_from_splitmappings.sh
chim2=$bashDir/find_chimeric_junctions_from_exon_to_exon_connections.sh
findGeneConnections=$bashDir/find_gene_to_gene_connections_from_pe_rnaseq.sh

# Bin 
gemrnatools=$binDir/gemtools-1.7.1-i3/bin/gem-rna-tools
bedtools=$binDir/bedtools2-2.20.1/bin/bedtools
gtfilter=$binDir/gemtools-1.7.1-i3/bin/gt.filter

# Awk 
unmapped=$awkDir/extract_unmapped.awk
bed12ToGff=$awkDir/bed12fields2gff.awk
gff2Gff=$awkDir/gff2gff.awk
gemToGff=$awkDir/gemsplit2gff_unique4.awk
bedCorrectStrand=$awkDir/bedCorrectStrand.awk
mapCorrectStrand=$awkDir/gemCorrectStrand.awk
addPEinfo=$awkDir/add_PE_info.awk
AddSimGnPairs=$awkDir/add_sim_bt_gnPairs.awk
juncFilter=$awkDir/chimjunc_filter.awk

## DISPLAY PIPELINE CONFIGURATION  
##################################


printf "\n"
header="Pipeline configuration for $lid"
echo $header
eval "for i in {1..${#header}};do printf \"-\";done"
printf "\n\n"
printf "  %-34s %s\n\n" "ChimPipe Version $version"
printf "  %-34s %s\n" "***** Mandatory *****"
printf "  %-34s %s\n" "Input file:" "$input"
printf "  %-34s %s\n" "Reference genome file:" "$index"
printf "  %-34s %s\n" "Reference gene annotation file:" "$annot"
printf "  %-34s %s\n" "Quality offset:" "$quality"
printf "  %-34s %s\n\n" "Sample identifier:" "$lid"

printf "  %-34s %s\n" "***** Reads information *****"
printf "  %-34s %s\n" "Strandedness (1:stranded,0:unstranded):" "$stranded"
printf "  %-34s %s\n" "Reads directionality:" "$readDirectionality"
printf "  %-34s %s\n\n" "Maximum read length:" "$maxReadLength"

printf "  %-34s %s\n" "***** Mapping *****"
printf "  %-34s %s\n" "Max number of allowed mismatches contiguous mapping:" "$mism"
printf "  %-34s %s\n" "Max number of allowed mismatches split mapping:" "$mismSplit"
printf "  %-34s %s\n" "Consensus-splice-sites:" "$spliceSites"
printf "  %-34s %s\n" "Minimum split size for segmental mapping:" "$splitSize"
printf "  %-34s %s\n\n" "Mapping statistics (1:enabled,0:disabled):" "$mapStats"

printf "  %-34s %s\n" "***** Filters *****"
printf "  %-34s %s\n" "Chimeric junctions filtering module configuration:" "$filterConf"
printf "  %-34s %s\n\n" "Similarity information between gene pairs:" "$simGnPairs"

printf "  %-34s %s\n" "***** General *****"
printf "  %-34s %s\n" "Output directory:" "$outDir"
printf "  %-34s %s\n" "Temporary directory:" "$TMPDIR"
printf "  %-34s %s\n" "Number of threads:" "$threads"
printf "  %-34s %s\n\n" "Loglevel:" "$logLevel"

## START CHIMPIPE
#################
header="Executing ChimPipe $version for $lid"
echo $header
eval "for i in {1..${#header}};do printf \"-\";done"
printf "\n\n"
pipelineStart=$(date +%s)

exit 1
# 1) map all the reads to the genome, to the transcriptome and de-novo, using the standard RNA-Seq 
#################################################################################################
#   mapping pipeline but with max intron size larger than the biggest chromosome, and with an edit
##################################################################################################  
#   distance of round(read_length/20) -> get a gem map file and a bam file of "normal" mappings;
################################################################################################
# outputs are: 
##############
# - $outDir/$lid.map.gz 
# - $outDir/$lid\_filtered_cuff.bam

bamFirstMapping=$outDir/${lid}_filtered_cuff.bam
if [ ! -e $bamFirstMapping ]; then
	step="FIRST-MAP"
	startTime=$(date +%s)
	printHeader "Executing first mapping step"    
    run "$pipeline -i $input -g $index -a $annot -q $quality -M $mism $mapStats -L $maxReadLength -e $lid -l $logLevel -t $threads" "$ECHO" 
	endTime=$(date +%s)
	printHeader "First mapping for $lid completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "First mapping BAM file is present...skipping first mapping step"
fi

# 2) extract the reads that do not map with a number of mismatches lower than 6
###############################################################################
#    from the gem file: get a fastq file 
########################################
# output is: 
############
# - $outDir/$lid.unmapped.fastq

unmappedReads=$outDir/${lid}.unmapped.fastq

if [ ! -e $unmappedReads ]; then
	step="UNMAP"
	startTime=$(date +%s)
	printHeader "Executing unmapped reads step" 
	log "Extracting the reads that do not map with a number of mismatches lower than 6..." $step
	run "zcat $lid.map.gz | awk -f $unmapped | $gtfilter -t $threads --output-format 'FASTA' > $outDir/$lid.unmapped.fastq 2> $outDir/$lid.unmapped.err" "$ECHO" 	
	log "done\n" 
    if [ -e $unmappedReads ]; then
    	log "Computing md5sum for unmapped reads file..." $step
    	run "md5sum $unmappedReads > $unmappedReads.md5" "$ECHO"
        log "done\n"
    else
        log "Error extracting the unmapped reads\n" "ERROR" 
        exit -1
    fi
    endTime=$(date +%s)
	printHeader "Extracting unmapped reads completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Unmapped reads file is present...skipping extracting unmapped reads step"
fi

# 3) map the unmapped reads with the rna mapper binary (!!!parameters to think!!!): get a gem map file
######################################################################################################
#   of "atypical" mappings
########################
# output is: 
############
# - $outDir/SecondMapping/$lid.unmapped_rna-mapped.map

gemSecondMapping=$outDir/SecondMapping/${lid}.unmapped_rna-mapped.map

if [ ! -e $gemSecondMapping ]; then
	step="SECOND-MAP"
	startTime=$(date +%s)
	printHeader "Executing second mapping step"
	log "Mapping the unmapped reads with the rna mapper..." $step	
	run "$gemrnatools split-mapper -I $index -i $outDir/$lid.unmapped.fastq -q 'offset-$quality' -o $outDir/SecondMapping/$lid.unmapped_rna-mapped -t 10 -T $threads -c $spliceSites --min-split-size $splitSize > $outDir/SecondMapping/$lid.gem-rna-mapper.out 2> $outDir/SecondMapping/$lid.gem-rna-mapper.err" "$ECHO"
	log "done\n" 
	if [ -e $gemSecondMapping ]; then
    	log "Computing md5sum for the gem file with the second mappings..." $step
    	run "md5sum $gemSecondMapping > $gemSecondMapping.md5" "$ECHO"
        log "done\n"
    else
        log "Error in the second mapping\n" "ERROR" 
        exit -1
    fi
    endTime=$(date +%s)
	printHeader "Unmapped reads mapping step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
	printHeader "Second mapping GEM file is present...skipping extracting second mapping step"
fi

# 4) extract the reads mapping both uniquely and in 2 blocks from the bam file of "normal" mappings and convert in gff.gz
#########################################################################################################################
# output is: 
############
# - $outDir/FromFirstBam/$lid\_filtered_cuff_2blocks.gff.gz

gffFromBam=$outDir/FromFirstBam/${lid}_filtered_cuff_2blocks.gff.gz

if [ ! -e $gffFromBam ]; then
	step="FIRST-CONVERT"
	startTime=$(date +%s)
	printHeader "Executing conversion of the bam into gff step"
	log "Generating a ".gff.gz" file from the normal mappings containing the reads split-mapping both uniquely and in 2 blocks..." $step
	$bedtools bamtobed -i $outDir/$lid\_filtered_cuff.bam -bed12 | awk '$10==2' | awk -v readDirectionality=$readDirectionality -f $bedCorrectStrand | awk -f $bed12ToGff | awk -f $gff2Gff | gzip > $outDir/FromFirstBam/$lid\_filtered_cuff_2blocks.gff.gz
	log "done\n"
	if [ -e $gffFromBam ]; then
    	log "Computing md5sum for the gff file from the ".bam" containing the reads mapping both uniquely and in 2 blocks..." $step
    	run "md5sum $gffFromBam > $gffFromBam.md5" "$ECHO"
        log "done\n"
    else
        log "Error Generating the gff file\n" "ERROR" 
        exit -1
    fi
	endTime=$(date +%s)
	printHeader "Step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
	printHeader "Gff from first mapping BAM file is present... skipping conversion step"
fi

# 5) extract the reads mapping both uniquely and in 2 blocks from the map file of "atypical" mappings and convert in gff.gz
#########################################################################################################################
# output is: 
############
# - $outDir/FromSecondMapping/${lid}.unmapped_rna-mapped.gff.gz

gffFromMap=$outDir/FromSecondMapping/${lid}.unmapped_rna-mapped.gff.gz

if [ ! -e $gffFromMap ]; then
	step="SECOND-CONVERT"
	startTime=$(date +%s)
	printHeader "Executing conversion of the gem into gff step"
	log "Generating a ".gff.gz" file from the atypical mappings containing the reads split-mapping both uniquely and in 2 blocks..." $step
	run "awk -v readDirectionality=$readDirectionality -f $mapCorrectStrand $outDir/SecondMapping/$lid.unmapped_rna-mapped.map | awk -v rev=1 -f $gemToGff | awk -f $gff2Gff | gzip > $outDir/FromSecondMapping/${lid}.unmapped_rna-mapped.gff.gz" "$ECHO"
	log "done\n" 
	if [ -e $gffFromMap ]; then
    	log "Computing md5sum for the gff file from the atypical mappings containing the reads mapping both uniquely and in 2 blocks..." $step
    	run "md5sum $gffFromMap > $gffFromMap.md5" "$ECHO"
        log "done\n"
    else
        log "Error Generating the gff file\n" "ERROR" 
        exit -1
    fi
    endTime=$(date +%s)
	printHeader "Step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
	printHeader "Gff from second mapping MAP file is present... skipping conversion step"
fi


# 6) put the path to the "normal" and "atypical" gff.gz files in a same txt file for chimsplice 
#############################################################################################

run "echo $outDir/FromFirstBam/$lid\_filtered_cuff_2blocks.gff.gz > $outDir/split_mapping_file_sample_$lid.txt" "$ECHO"
run "echo $outDir/FromSecondMapping/$lid.unmapped_rna-mapped.gff.gz >> $outDir/split_mapping_file_sample_$lid.txt" "$ECHO"

# 7) run chimsplice on 4) and 5)
###############################
# - $outDir/Chimsplice/chimeric_junctions_report_$lid.txt
# - $outDir/Chimsplice/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt
# - $outDir/Chimsplice/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn_morethan10staggered.txt

exonConnections1=$outDir/Chimsplice/exonA_exonB_with_splitmapping_part1overA_part2overB_readlist_sm1list_sm2list_staggeredlist_totalist_$lid\_filtered_cuff_2blocks.gff.txt.gz 
exonConnections2=$outDir/Chimsplice/exonA_exonB_with_splitmapping_part1overA_part2overB_readlist_sm1list_sm2list_staggeredlist_totalist_$lid.unmapped_rna-mapped.gff.txt.gz

printHeader "Executing Chimsplice step"
if [ ! -e $exonConnections1 ] || [ ! -e $exonConnections2 ]; then
	step="CHIMSPLICE"
	startTime=$(date +%s)
	log "Finding exon to exon connections from the ".gff.gz" files containing the "normal" and "atypical" mappings..." $step
	run "$chim1 $outDir/split_mapping_file_sample_$lid.txt $annot $outDir/Chimsplice $stranded 2> $outDir/Chimsplice/find_exon_exon_connections_from_splitmappings_better_$lid.err" "$ECHO"
	log "done\n" 
	if [ ! -e $exonConnections1 ] || [ ! -e $exonConnections2 ]; then
        log "Error running chimsplice\n" "ERROR" 
        exit -1
    fi
    endTime=$(date +%s)
	printHeader "Find exon to exon connections step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
	printHeader "Exon to exon connections file present... skipping step"
fi

chimJunctions=$outDir/Chimsplice/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt
if [ ! -e $chimJunctions ]; then
	step="CHIMSPLICE"
	startTime=$(date +%s)
	log "Finding chimeric junctions from exon to exon connections..." $step
	run "$chim2 $outDir/split_mapping_file_sample_$lid.txt $index $annot $outDir/Chimsplice $stranded > $outDir/Chimsplice/chimeric_junctions_report_$lid.txt 2> $outDir/Chimsplice/find_chimeric_junctions_from_exon_to_exon_connections_$lid.err" "$ECHO"
	log "done\n" 
	if [ ! -e $chimJunctions ]; then
        log "Error running chimsplice\n" "ERROR" 
        exit -1
    fi
    endTime=$(date +%s)
	printHeader "Find chimeric junctions step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
	printHeader "Chimeric Junctions file present... skipping step"
fi

# 8.1) Find gene to gene connections supported by paired-end mappings from the bam file of "normal" mappings with the number of mappings supporting the connection.  
#################################################################################################################################################################
# For a connection g1 to g2 to exist there must be at least one mapping where the first mate is strandedly (if data is stranded) overlapping with an exon 
#########################################################################################################################################################
# of g1 and the second mate is (strandedly if data is stranded) overlapping with an exon of g2
##############################################################################################
# - $outDir/PE/readid_gnlist_whoseexoverread_noredund.txt.gz
# - $outDir/PE/readid_twomateswithgnlist_alldiffgnpairs_where_1stassociatedto1stmate_and2ndto2ndmate.txt.gz
# - $outDir/PE/pairs_of_diff_gn_supported_by_pereads_nbpereads.txt

PEinfo=$outDir/PE/pairs_of_diff_gn_supported_by_pereads_nbpereads.txt

printHeader "Executing find gene to gene connections from PE mappings step"
if [ ! -e $PEinfo ];then
	step="PAIRED-END"
	startTime=$(date +%s)
	log "Finding gene to gene connections supported by paired-end mappings from the ".bam" containing reads mapping in a unique and continuous way..." $step
	run "$findGeneConnections $outDir/$lid\_filtered_cuff.bam $annot $outDir/PE $readDirectionality 2> $outDir/PE/find_gene_to_gene_connections_from_pe_rnaseq_fast_$lid.err" "$ECHO"
	log "done\n" 
	if [ ! -e $PEinfo ]; then
        log "Error finding gene to gene connections\n" "ERROR" 
        exit -1
    fi
    endTime=$(date +%s)
	printHeader "Find gene to gene connections from PE mappings step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
	printHeader "Gene to gene connections file present... skipping step"
fi

# 8.2) Add gene to gene connections information to chimeric junctions matrix
##########################################################################
# - $outDir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_PEinfo_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt

chimJunctionsPE=$outDir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_PEinfo_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt

if [ ! -e $chimJunctionsPE ];then
	step="PAIRED-END"
	startTime=$(date +%s)
	log "Adding PE information to the matrix containing chimeric junction candidates..." $step
	run "awk -v fileRef=$PEinfo -f $addPEinfo $chimJunctions 1> $outDir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_PEinfo_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt
" "$ECHO"
	log "done\n" 
	if [ ! -e $chimJunctionsPE ]; then
        log "Error adding PE information\n" "ERROR" 
        exit -1
    fi
    endTime=$(date +%s)
	printHeader "Add PE information step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
	printHeader "Chimeric junction matrix with PE information present... skipping step"
fi

# 9) Add information regarding the sequence similarity between connected genes
##############################################################################
# - $outDir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_PEinfo_maxLgalSim_maxLgal_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt

chimJunctionsSim=$outDir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_PEinfo_maxLgalSim_maxLgal_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt

if [  ! -e "$chimJunctionsSim" ]; then		
	if [ "$simGnPairs" != "" ];then
		step="SIM"
		startTime=$(date +%s)
		log "Adding sequence similarity between connected genes information to the chimeric junction matrix..." $step
		run "awk -v fileRef=$simGnPairs -f $AddSimGnPairs $chimJunctionsPE 1> $chimJunctionsSim" "$ECHO"
		log "done\n" 
		if [ ! -e $chimJunctionsSim ]; then
			log "Error adding similarity information\n" "ERROR" 
	    	exit -1
		fi
		rm $chimJunctionsPE
		endTime=$(date +%s)
		printHeader "Add sequence similarity information step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
	else 
		printHeader "Similarity information between the gene pairs in the annotation does not provided... skipping step"
	fi
else
	printHeader "Chimeric junction matrix with similarity information present... skipping step"
fi

# 10) Produce a matrix containing chimeric junction candidates with a header in the first row
#############################################################################################
# - $outDir/chimeric_junctions_candidates.txt

chimJunctionsCandidates=$outDir/chimeric_junctions_candidates.txt

if [ ! -e "$chimJunctionsCandidates" ]; then		
	step="HEADER"
	log "Adding a header to the matrix containing the chimeric junction candidates..." $step
	if [ -e "$chimJunctionsSim" ]; then		
		run "awk 'BEGIN{print "juncId", "nbstag", "nbtotal", "maxbeg", "maxEnd", "samechr", "samestr", "dist", "ss1", "ss2", "gnlist1", "gnlist2", "gnname1", "gnname2", "bt1", "bt2", "PEsupport", "maxSim", "maxLgal";}{print \$0;}' $chimJunctionsSim 1> $chimJunctionsCandidates" "$ECHO"
		rm $chimJunctionsSim
		log "done\n"
	else
		if [ -e "$chimJunctionsPE" ]; then
			run "awk 'BEGIN{print "juncId", "nbstag", "nbtotal", "maxbeg", "maxEnd", "samechr", "samestr", "dist", "ss1", "ss2", "gnlist1", "gnlist2", "gnname1", "gnname2", "bt1", "bt2", "PEsupport";}{print \$0;}' $chimJunctionsPE 1> $chimJunctionsCandidates" "$ECHO"
			log "done\n" 
			rm $chimJunctionsPE
		else
			log "Error, intermediate file: $chimJunctionsSim or $chimJunctionsPE is missing\n" "ERROR" 
	    	exit -1			
		fi
	fi 
else
	printHeader "Header already added... skipping step"
fi

# 11) Filter out chimera candidates to produce a final set of chimeric junctions
################################################################################
# - $outDir/chimeric_junctions.txt

chimJunctions=$outDir/chimeric_junctions.txt

if [ ! -e "$chimJunctions" ]; then
	step="FILTERING MODULE"
	startTime=$(date +%s)
	log "Filtering out chimera candidates to produce a final set of chimeric junctions..." $step
	awk -v filterConf=$filterConf -f $juncFilter $chimJunctionsCandidates 1> $chimJunctions
	log "done\n" 
	if [ ! -e $chimJunctionsSim ]; then
		log "Error filtering chimeric junction candidates\n" "ERROR" 
    	exit -1
	fi
	rm $chimJunctionsCandidates
	endTime=$(date +%s)
	printHeader "Filtering module step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else 
	printHeader "Chimeric junction candidates already filtered... skipping step"
fi

# 11) END
#########
pipelineEnd=$(date +%s)
printHeader "Chimera Mapping pipeline for $lid completed in $(echo "($pipelineEnd-$pipelineStart)/60" | bc -l | xargs printf "%.2f\n") min "
exit 0

