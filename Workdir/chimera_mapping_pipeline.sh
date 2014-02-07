#!/bin/bash
#
#Forcing bash shell
#
#$ -S /bin/bash
#
#$ -V
#
#$ -cwd
# -M rodriguezmartinbernardo@gmail.com
# -m b
#
#$ -pe smp 8
#$ -q long,rg-el6
#$ -l virtual_free=64G
#$ -o $JOB_NAME.out
#$ -e $JOB_NAME.err

# will exit if there is an error or in a pipe
set -e -o pipefail

function usage
{
cat <<instructions
USAGE: chimera_mapping_pipeline.sh -f <fastq_file> -i <genome_index> -a <annotation> -q <quality> [options]
       
Execute Chimera mapping pipeline (from fastq file to chimeric junctions detection) on one sample

** Mandatory arguments:

	-f	<INPUT_FILE>	First mate FASTQ file. Pipeline designed to deal with paired end data. 
				Please make sure the second mate is in the same directory and the files are named according to the following: 
				"YourSampleId"_1.fastq.gz and "YourSampleId"_2.fastq.gz; where "YourSampleId" must be the same for the 
				mate one (_1) than for the mate two (_2)  
	
	-i	<GEM>		Index for the reference genome (".gem" format).
	-a	<GTF>		Reference annotation (".gtf" format).
	-q	<NUMBER>	Quality offset of the reads in the ".fastq" [33 | 64 | ignore].
	-e	<STRING> 	Sample identifier (the output files will be named according this id).

** [options] can be:
 
	-b			Flag to specify that the input file has BAM format (Already mapped reads). Not enabled 
	-s			Flag to specify whether the sample has strandness information for the reads. Default false (So data unstranded).
	-d	<STRING>	Directionality of the reads (MATE1_SENSE, MATE2_SENSE, MATE_STRAND_CSHL, SENSE, ANTISENSE & NONE). Default "NONE".
	-m	<NUMBER>	Max number of mismatches. Default 4.
	-M	<NUMBER>	Max read length. This is used to create the de-novo transcriptome and acts as an upper bound. Default 150.
	-c	<pair_1>, ... ,<pair_s>
      		with <pair> := <donor_consensus>+<acceptor_consensus> (list of pairs of donor/acceptor splice site consensus sequences. Default "GT+AG")
                          	 
	-S	<NUMBER>	Minimum split size. Default 15
	-o	<PATH>		Output directory (By default it uses the current working directory).
	
	-l	<STRING>	Log level (error, warn, info, debug). Default "info".
	-h			Flag to display usage information.
	-t			Flag to test the pipeline. Writes the commands to the standard output.
	exit 0
instructions
}


function printHeader {
    string=$1
    echo "`date` ***** $string *****"
}

function log {
    string=$1
    label=$2
    if [[ ! $ECHO ]];then
        if [[ $label != "" ]];then
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


# PARSING INPUT ARGUMENTS 
#########################
while getopts ":f:i:a:q:e:bsd:m:M:c:S:o:l:th" opt; do
  case $opt in
    f)
      input="$OPTARG"
      ;;
    i)
      index="`dirname $OPTARG | xargs readlink -e`/`basename $OPTARG`"
      ;;
    a)
      annot="`dirname $OPTARG | xargs readlink -e`/`basename $OPTARG`"
      ;;
    q)
      quality=$OPTARG
      ;;
    e)
      lid=$OPTARG
      ;;     
    b)
      bam=1
      ;;
    s)
      stranded=1
      ;;
    d)
      readDirectionality=$OPTARG
      ;;
    m)
      mism=$OPTARG
      ;;
    M)
      maxReadLength=$OPTARG 
      ;; 
 	c)
 	  spliceSites=$OPTARG
 	  ;;
 	S)
 	  splitSize=$OPTARG
 	  ;;  
 	o)
      outDir=$OPTARG
      ;;
 	l)
      logLevel=$OPTARG
      ;;
    t)
      ECHO="echo "
      ;;
    h)
      usage
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# SETTING VARIABLES AND INPUT FILES
###################################
if [[ ! -e $input ]]; then log "Please specify a valid input file\n" "ERROR" >&2; exit -1; fi
if [[ `basename ${input##*_}` != "1.fastq.gz" ]]; then log "Please check that the name of your FASTQ file ends with \"_1.fastq.gz\"\n" "ERROR" >&2; exit -1; fi
if [[ ! -e $index ]]; then log "Please specify a valid genome index file\n" "ERROR" >&2; exit -1; fi
if [[ ! -e $annot ]]; then log "Please specify a valid annotation file\n" "ERROR" >&2; exit -1; fi
if [[ $quality == "" ]]; then log "Please specify the quality\n" "ERROR" >&2; exit -1; fi
if [[ $lid == "" ]]; then log "Please specify the experiment identifier\n" "ERROR" >&2; exit -1; fi
if [[ $bam == "" ]]; then bam=0; fi
if [[ $stranded == "" ]]; then stranded=0; fi
if [[ $readDirectionality == "" ]]; then readDirectionality='NONE'; fi
if [[ $mism == "" ]]; then mism='4'; fi
if [[ $maxReadLength == "" ]]; then maxReadLength='150'; fi
if [[ $spliceSites ==  "" ]]; then spliceSites='GT+AG' ; fi
if [[ $splitSize == "" ]]; then splitSize='15'; fi
if [[ ! -d $outDir ]]; then outDir=${SGE_O_WORKDIR-$PWD}; fi
if [[ $logLevel == "" ]]; then logLevel='info'; fi


# SETTING UP THE ENVIRONMENT
############################

# = Directories = #
binDir=~brodriguez/Chimeras_project/Chimeras_detection_pipeline/Chimera_mapping/Workdir/bin
awkDir=~brodriguez/Chimeras_project/Chimeras_detection_pipeline/Chimera_mapping/Workdir/Awk
bashDir=~brodriguez/Chimeras_project/Chimeras_detection_pipeline/Chimera_mapping/Workdir/Bash
if [[ ! -d $outDir/SecondMapping ]]; then mkdir $outDir/SecondMapping; fi
if [[ ! -d $outDir/FromFirstBam ]]; then mkdir $outDir/FromFirstBam; fi
if [[ ! -d $outDir/FromSecondMapping ]]; then mkdir $outDir/FromSecondMapping; fi
if [[ ! -d $outDir/Chimsplice ]]; then mkdir $outDir/Chimsplice; fi

# = Programs/Scripts = #
# Bash 
pipeline=~brodriguez/Chimeras_project/Chimeras_detection_pipeline/Chimera_mapping/Versions/V0.3.2/blueprint.pipeline.sh 
chim1=~brodriguez/Chimeras_project/Chimeras_detection_pipeline/Chimsplice/Versions/V0.3.0/find_exon_exon_connections_from_splitmappings_better2.sh
chim2=~brodriguez/Chimeras_project/Chimeras_detection_pipeline/Chimsplice/Versions/V0.3.0/find_chimeric_junctions_from_exon_to_exon_connections_better2.sh

# Bin 
rnaMapper=$binDir/gem-rna-mapper
bamToBed=$binDir/bamToBed

# Awk 
bed12ToGff=$awkDir/bed12fields2gff.awk
gff2Gff=$awkDir/gff2gff.awk
gemToGff=$awkDir/gemsplit2gff_unique4.awk
bedCorrectStrand=$awkDir/bedCorrectStrand.awk
mapCorrectStrand=$awkDir/gemCorrectStrand.awk

# Python 
unmapped=$binDir/filter_unmapped.py 

# Activating gemtools environment
source ~sdjebali/ENCODE_AWG/Analyses/Mouse_Human/Chimeras/gemtools/environment/bin/activate


## DISPLAY PIPELINE CONFIGURATION  
##################################

printf "\n"
printf "*****Chimera Mapping pipeline configuration*****\n"
printf "Pipeline Version: V0.1.1\n"
printf "Input: $input\n"
printf "Index: $index\n"
printf "Annotation: $annot\n"
printf "Quality: $quality\n"
printf "Strand: $stranded\n"
printf "Directionality: $readDirectionality\n"
printf "Number mismatches: $mism\n"
printf "Max read length: $maxReadLength\n"
printf "Consensus split sites: $spliceSites\n"
printf "Min split size: $splitSize\n"
printf "Outdir: $outDir\n"
printf "Exp. id: $lid\n"
printf "Log level: $logLevel\n"	
printf "\n"

## START CHIMERA MAPPING PIPELINE 
#################################

printHeader "Starting Chimera Mapping pipeline for $lid"
pipelineStart=$(date +%s)

# 1) map all the reads to the genome, to the transcriptome and de-novo, using the standard rnaseq 
#################################################################################################
#   mapping pipeline but with max intron size larger than the biggest chromosome, and with an edit
##################################################################################################  
#   distance of round(read_length/20) -> get a gem map file and a bam file of "normal" mappings;
################################################################################################
# outputs are: 
##############
# - $outDir/$lid.map.gz 
# - $outDir/$lid\_filtered_cuff.bam
printHeader "Executing first mapping step with the Blueprint pipeline"

run "$pipeline -i $input -g $index -a $annot -q $quality -m $mism -M $maxReadLength -e $lid -l $logLevel" "$ECHO"

# 2) extract the reads that do not map with an edit distance strictly greater than round(read_length/20) 
####################################################################################################
#    from the gem file: get a fastq file
########################################
# output is: 
############
# - $outDir/$lid.unmapped.fastq
startTime=$(date +%s)
printHeader "Extracting the reads that do not map with an edit distance strictly greater than round (read_length/20)"

run "$unmapped -i $lid.map.gz -t 8 -f fastq -m 5 > $lid.unmapped.fastq 2> $lid.unmapped.err" "$ECHO"

endTime=$(date +%s)
printHeader "Extracting unmapped reads completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"

# 3) map the unmapped reads with the rna mapper binary (!!!parameters to think!!!): get a gem map file
######################################################################################################
#   of "exotic" mappings
########################
# output is: 
############
# - $outDir/SecondMapping/$lid.unmapped_rna-mapped.map
startTime=$(date +%s)
printHeader "Mapping the unmapped reads with the rna mapper"

run "$rnaMapper -I $index -i $outDir/$lid.unmapped.fastq -q 'offset-$quality' -o $outDir/SecondMapping/$lid.unmapped_rna-mapped -t 10 -T 8 -c $spliceSites --min-split-size $splitSize > $outDir/SecondMapping/$lid.gem-rna-mapper.out 2> $outDir/SecondMapping/$lid.gem-rna-mapper.err" "$ECHO"

endTime=$(date +%s)
printHeader "Unmapped reads mapping step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"

# 4) extract the reads mapping both uniquely and in 2 blocks from the bam file of "normal" mappings and convert in gff.gz
#########################################################################################################################

startTime=$(date +%s)
printHeader "Generating a ".gff.gz" file from the ".bam" containing the reads mapping both uniquely and in 2 blocks"

$bamToBed -i $outDir/$lid\_filtered_cuff.bam -bed12 | awk '$10==2' | awk -v readDirectionality=$readDirectionality -f $bedCorrectStrand | awk -f $bed12ToGff | awk -f $gff2Gff | gzip > $outDir/FromFirstBam/$lid\_filtered_cuff_2blocks.gff.gz

endTime=$(date +%s)
printHeader "Step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"

# 5) extract the reads mapping both uniquely and in 2 blocks from the map file of "exotic" mappings and convert in gff.gz
#########################################################################################################################

startTime=$(date +%s)
printHeader "Generating a ".gff.gz" file from the ".bam" containing the "exotic" mappings"

run "awk -v readDirectionality=$readDirectionality -f $mapCorrectStrand $outDir/SecondMapping/$lid.unmapped_rna-mapped.map | awk -v rev=1 -f $gemToGff | awk -f $gff2Gff | gzip > $outDir/FromSecondMapping/$lid.unmapped_rna-mapped.gff.gz" "$ECHO"

endTime=$(date +%s)
printHeader "Step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"

# 6) put the path to the "normal" and "exotic" gff.gz files in a same txt file for chimsplice 
#############################################################################################

run "echo $outDir/FromFirstBam/$lid\_filtered_cuff_2blocks.gff.gz > $outDir/split_mapping_file_exp_$lid.txt" "$ECHO"
run "echo $outDir/FromSecondMapping/$lid.unmapped_rna-mapped.gff.gz >> $outDir/split_mapping_file_exp_$lid.txt" "$ECHO"

# 7) run chimsplice on 4) and 5)
###############################
startTime=$(date +%s)
printHeader "Running Chimsplice on the ".gff.gz" files containing the "normal" and "exotic" mappings "

run "$chim1 $outDir/split_mapping_file_exp_$lid.txt $annot $outDir/Chimsplice $stranded 2> $outDir/Chimsplice/find_exon_exon_connections_from_splitmappings_better_$lid.err" "$ECHO"

run "$chim2 $outDir/split_mapping_file_exp_$lid.txt $annot $outDir/Chimsplice $stranded > $outDir/Chimsplice/chimeric_junctions_report_$lid.txt 2> $outDir/Chimsplice/find_chimeric_junctions_from_exon_to_exon_connections_$lid.err" "$ECHO"

endTime=$(date +%s)
printHeader "Chimsplice step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"


# 8) END
########

pipelineEnd=$(date +%s)

printHeader "Chimera Mapping pipeline for $lid completed in $(echo "($pipelineEnd-$pipelineStart)/60" | bc -l | xargs printf "%.2f\n") min "
exit 0

