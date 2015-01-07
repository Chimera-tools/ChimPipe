#!/bin/bash

<<authors
*****************************************************************************
	
	ChimPipe.sh
	
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

# FUNCTIONS
############

# Function 1. Print basic usage information
############################################
function usage
{
cat <<help
	
**** ChimPipe version $version ****

Execute ChimPipe on one Illumina paired-end RNA-seq dataset (sample).
	
*** USAGE

FASTQ:

	$0 --fastq_1 <mate1_fastq> --fastq_2 <mate2_fastq> -g <genome_index> -a <annotation> -t <transcriptome_index> -k <transcriptome_keys> [OPTIONS]

BAM:	

	$0 --bam <bam> -g <genome_index> -a <annotation> [OPTIONS]

*** MANDATORY ARGUMENTS
		
* FASTQ:

	--fastq_1			<FASTQ>		First mate sequencing reads in FASTQ format. It can be gzip compressed [.gz].
	--fastq_2			<FASTQ>		Second mate sequencing reads in FASTQ format. It can be gzip compressed [.gz].
	-g|--genome-index		<GEM>		Reference genome index in GEM format.
	-a|--annotation			<GTF>		Reference gene annotation file in GTF format.                                			
	-t|--transcriptome-index	<GEM>		Annotated transcriptome index in GEM format.
	-k|--transcriptome-keys		<KEYS>		Transcriptome to genome coordinate conversion keys.  
	--sample-id			<STRING>	Sample identifier (output files are named according to this id).  
	
* BAM:	

	--bam				<BAM>		Mapped reads in BAM format. A splicing aware aligner is needed to map the reads. 
	-g|--genome-index		<GEM>		Reference genome index in GEM format.
	-a|--annotation			<GTF>		Reference genome annotation file in GTF format.
	--sample-id			<STRING>	Sample identifier (the output files will be named according to this id).  
	
*** [OPTIONS] can be:

* General: 
	--log				<STRING>	Log level [error | warn | info | debug]. Default warn.
	--threads			<INTEGER>	Number of threads to use. Default 1.
	-o|--output-dir			<PATH>		Output directory. Default current working directory. 
	--tmp-dir			<PATH>		Temporary directory. Default /tmp.	
	--no-cleanup					Keep intermediate files. 		
	--dry						Test the pipeline. Writes the commands to the standard output.
	-h|--help					Display partial usage information, only mandatory plus general arguments.
	-f|--full-help					Display full usage information with additional options. 

help
}


# Function 2. Print all the other options
##########################################
function usage_long
{
cat <<help
* Read information:
	--max-read-length		<INTEGER>	Maximum read length. This is used to create the de-novo transcriptome and acts as an upper bound. Default 150.
	-l|--seq-library 		<STRING> 	Type of sequencing library [MATE1_SENSE | MATE2_SENSE | UNSTRANDED].
                        				UNSTRANDED for not strand-specific protocol (unstranded data) and the others for the different types 
							of strand-specific protocols (stranded data).
* Mapping phase parameters
                        				
  First mapping:
	-C|--consensus-ss-fm		<(couple_1)>, ... ,<(couple_s)>	with <couple> := <donor_consensus>+<acceptor_consensus>
                                 			List of couples of donor/acceptor splice site consensus sequences. Default='GT+AG,GC+AG,ATATC+A.,GTATC+AT'
	-S|--min-split-size-fm		<INTEGER>	Minimum split size for the segmental mapping steps. Default 15.
	--refinement-step-size-fm   	<INTEGER>   	If not mappings are found a second attempt is made by eroding "N" bases toward the ends of the read. 
							A value of 0 disables it. Default 2. 
	--no-stats					Disable mapping statistics. Default enabled.

  Second Mapping:
	--no-remap-unmapped				No remap first mapping unmapped reads. Default remapped 
	--remap-multimapped		<ALL|INTEGER> 	Remap first mapping multimapped reads (ALL: remap all the multimapped reads; 
							INTEGER: remap multimapped reads with more than N matches). 
							Default multimapped not remapped. 
	--remap-unique 			<INTEGER> 	Remap first mapping uniquelly mapped reads with more than N mismatches (Mismatches 
							in bad-quality bases also considered). Default unique mappings not remapped
	-c|--consensus-ss-sm		<(couple_1)>, ... ,<(couple_s)>	List of couples of donor/acceptor splice site consensus sequences. Default='GT+AG'
	-s|--min-split-size-sm		<INTEGER>	Minimum split size for the segmental mapping steps. Default 15.
	--refinement-step-size-sm   	<INTEGER>   	If not mappings are found a second attempt is made by eroding "N" bases toward the ends of the read. 
							A value of 0 disables it. Default 2. 
    
* Chimera detection phase parameters:

	--consider-multimapped 		<ALL|INTEGER>	Consider multimapped reads for the detection of chimeric junctions (ALL: all the multimapped reads; 
							INTEGER: multimapped reads with N or less matches). Default multimapped not considered. 
	--similarity-gene-pairs	<TEXT>			Text file with similarity information between the gene pairs in the annotation.
							Needed for the filtering module to discard chimeric junctions connecting highly similar genes. 
							If this file is not provided it will be computed by ChimPipe.
													
help
}

# Function 3. Display a link to ChimPipe's documentation
########################################################
function doc
{
cat <<help
A complete documentation for ChimPipe can be found at: http://chimpipe.readthedocs.org/en/latest/index.html		
help
}

# Function 4. Short help 
#########################
function usagedoc
{
usage
doc
}

# Function 5. Long help 
#########################
function usagelongdoc
{
usage
usage_long
doc
}

# Function 6. Print a section header for the string variable
##############################################################
function printHeader {
    string=$1
    echo "`date` ***** $string *****"
}

# Function 7. Print a subsection header for the string variable
################################################################
function printSubHeader {
    string=$1
    echo "`date` * $string *"
}

# Function 8. Print log information (Steps and errors)
#######################################################
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

# Function 6. Execute and print to stdout commands 
###################################################
function run {
    command=($1)
    if [[ $2 ]];then
         ${2}${command[@]}
    else
        echo -e "\n"" "${command[@]}
        eval ${command[@]}
    fi
}

# Function 7. Copy annotation, genome index, transcriptome index 
#################################################################
# and keys to the temporary directory
#####################################
function copyToTmp {
    IFS=',' read -ra files <<< "$1"
    for i in ${files[@]};do
        case $i in
            "annotation")
                if [[ ! -e $TMPDIR/`basename $annot` ]];then
                    log "Copying annotation file to $TMPDIR..." $step
                    run "cp $annot $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "index")
                if [[ ! -e $TMPDIR/`basename $genomeIndex` ]];then
                    log "Copying genome index file to $TMPDIR..." $step
                    run "cp $genomeIndex $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "t-index")
                if [[ ! -e $TMPDIR/`basename $transcriptomeIndex` ]];then
                    log "Copying annotated transcriptome index file to $TMPDIR..." $step
                    run "cp $transcriptomeIndex $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "keys")
                if [[ ! -e $TMPDIR/`basename $transcriptomeKeys` ]];then
                    log "Copying annotated transcriptome keys file to $TMPDIR..." $step
                    run "cp $transcriptomeKeys $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            esac
    done
}

# Function 8. Run the gemtools RNA-Seq pipeline to map all the reads  
###################################################################
# to the genome, to the transcriptome and de-novo
##################################################
# Input files:
# - $fastq1 
# - $fastq2
# - $genomeIndex
# - $annot
# - $transcriptomeIndex
# - $transcriptomeKeys
# Output files:
# -	${lid}_firstMap.map.gz
# - ${lid}_firstMap.bam
# - ${lid}_firstMap.bam.bai
# - ${lid}_firstMap.stats.txt
# - ${lid}_firstMap.stats.json


function runGemtoolsRnaPipeline {

	gemFirstMap=$outDir/${lid}_firstMap.map.gz
	baiFirstMap=$outDir/${lid}_firstMap.bam.bai
	statsFirstMap=$outDir/${lid}_firstMap.stats.txt
	statsJsonFirstMap=$outDir/${lid}_firstMap.stats.json
	
	step="FIRST-MAP"
	startTimeFirstMap=$(date +%s)
	printHeader "Executing first mapping step"    
	    
	## Copy needed files to TMPDIR
    copyToTmp "index,annotation,t-index,keys"

    log "Running gemtools rna pipeline on ${lid}..." $step
    run "$gemtools --loglevel $logLevel rna-pipeline -f $fastq1 $fastq2 -i $TMPDIR/`basename $genomeIndex` -a $TMPDIR/`basename $annot` -r $TMPDIR/`basename $transcriptomeIndex` -k $TMPDIR/`basename $transcriptomeKeys` -q $quality --max-read-length $maxReadLength --max-intron-length 300000000 --min-split-size $splitSizeFM --refinement-step $refinementFM --junction-consensus $spliceSitesFM --no-filtered --no-xs $stats --no-count -n `basename ${bamFirstMap%.bam}` --compress-all --output-dir $TMPDIR -t $threads" "$ECHO" 
	log "done\n"
   
   	if [ -s $TMPDIR/`basename $bamFirstMap` ]; 
   	then 
	    # Checksums
	    log "Computing md5sum for map file..." $step
       	run "md5sum $TMPDIR/`basename $gemFirstMap` > ${gemFirstMap}.md5" "$ECHO"
       	log "Computing md5sum for bam file..." $step
       	run "md5sum $TMPDIR/`basename $bamFirstMap` > ${bamFirstMap}.md5" "$ECHO"
		
		# Copy files from temporary to output directory
		run "cp $TMPDIR/`basename $gemFirstMap` $gemFirstMap" "$ECHO"
       	run "cp $TMPDIR/`basename $bamFirstMap` $bamFirstMap" "$ECHO"
       	run "cp $TMPDIR/`basename $baiFirstMap` $baiFirstMap" "$ECHO"
       	
       	if [[ "$mapStats" == "true" ]]; 
		then 
       		run "cp $TMPDIR/`basename $statsFirstMap` $statsFirstMap" "$ECHO"
       		run "cp $TMPDIR/`basename $statsJsonFirstMap` $statsJsonFirstMap" "$ECHO"
   		fi
   	else
       	log "Error producing bam file\n" "ERROR"
       	exit -1
   	fi
   	endTimeFirstMap=$(date +%s)        		
	printHeader "First mapping for $lid completed in $(echo "($endTimeFirstMap-$startTimeFirstMap)/60" | bc -l | xargs printf "%.2f\n") min"
}

# Function 9. Parse user's input
################################
function getoptions {
ARGS=`$getopt -o "g:a:t:k:o:hfl:C:S:c:s:" -l "fastq_1:,fastq_2:,bam:,genome-index:,annotation:,transcriptome-index:,transcriptome-keys:,sample-id:,log:,threads:,output-dir:,tmp-dir:,no-cleanup,dry,help,full-help,max-read-length:,seq-library:,consensus-ss-fm:,min-split-size-fm:,refinement-step-size-fm:,no-stats,no-remap-unmapped,remap-multimapped:,remap-unique:,consensus-ss-sm:,min-split-size-sm:,refinement-step-size-sm:,consider-multimapped:,filter-chimeras:,similarity-gene-pairs:" \
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
   	
   	## MANDATORY ARGUMENTS
    --fastq_1)
      if [ -n "$2" ];
      then
        fastq1=$2
      fi
      shift 2;;

	--fastq_2)
      if [ -n "$2" ];
      then
        fastq2=$2
      fi
      shift 2;;
      
    --bam)
      if [ -n "$2" ];
      then
        bam=$2
        bamAsInput="true"
      fi
      shift 2;;
        
    -g|--genome-index)
      if [ -n "$2" ];
      then
        genomeIndex=$2
      fi
      shift 2;;

    -a|--annotation)
      if [ -n "$2" ];
      then
        annot=$2
      fi
      shift 2;;
          
    -t|--transcriptome-index)
      if [ -n "$2" ];
      then
        transcriptomeIndex=$2
      fi
      shift 2;;
    
    -k|--transcriptome-keys)
      if [ -n "$2" ];
      then
        transcriptomeKeys=$2
      fi
      shift 2;;    
	
	--sample-id)
       if [ -n "$2" ];
       then
         lid=$2
       fi
       shift 2;;
       
    ## OPTIONS
    
	# General:
    
    --log)
       if [ -n $2 ];
       then
         logLevel=$2
       fi
       shift 2;;
	
    --threads)
       if [ -n $2 ];
       then
      	 threads=$2
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
 	
	--no-cleanup)
	    cleanup="false";
	    shift;;   	
 
	--dry)
	    ECHO="echo "
	    shift;;

	-h|--help)
	    usagedoc;
	    exit 1
	    shift;;
    
	-f|--full-help)
	    usagelongdoc;
	    exit 1
	    shift;;	
    
    # Reads information:
    --max-read-length)
      if [ -n $2 ];
      then
        maxReadLength=$2
      fi
      shift 2;;
      
    -l|--seq-library)
      if [ -n $2 ];
      then
        readDirectionality=$2
      fi
      shift 2;;
     
    # First mapping parameters: 
    -C|--consensus-ss-fm)
	    if [ -n "$2" ];
	    then
			spliceSitesFM=$2
	    fi
	    shift 2;;
   
    -S|--min-split-size-fm)
    	if [ -n "$2" ];
	    then
		    splitSizeFM=$2
	    fi
	    shift 2;;
	
	--refinement-step-size-fm)
	  if [ -n "$2" ];
	  then
	      refinementFM=$2
	  fi
	  shift 2;;
      
    --no-stats)
    	if [ -n "$2" ];
	 	then
	  		mapStats="false"  
	  	fi
	  	shift;;
      
	# Second mapping parameters:
	--no-remap-unmapped)
		if [ -n "$2" ];
	  	then
	    	remapUnmapped="false"
	  	fi
	  	shift;;
	 
	--remap-multimapped)
		if [ -n "$2" ];
	  	then
	    	remapMultimapped="true"
	    	nbMultimaps2Remap=$2
	  	fi
	  	shift 2;;
	
	--remap-unique)
	    if [ -n "$2" ];
	  	then
	    	remapUnique="true"
	    	nbMism=$2
	  	fi
	  	shift 2;;
	  		
	-c|--consensus-ss-sm)
	    if [ -n "$2" ];
	    then
			spliceSitesSM=$2
	    fi
	    shift 2;;
	
	-s|--min-split-size-sm)
    	if [ -n "$2" ];
	    then
	    splitSizeSM=$2
	fi
	    shift 2;;
	
	--refinement-step-size-sm)
		if [ -n "$2" ];
	    then
		    refinementSM=$2
	    fi
	    shift 2;;
	
	# Chimera detection phase parameters:
	--consider-multimapped)
		if [ -n "$2" ];
	  	then
	    	multimapped2ChimDetection="true"
	  		nbMultimaps2ChimDetect=$2
	  	fi
	  	shift 2;;
	  	
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
      
	--)
	    shift
	    break;;  
    esac
done
}

# SETTING UP THE ENVIRONMENT
############################

# ChimPipe version 
version=v0.8.8

# Enable extended pattern matching 
shopt -s extglob

# 1. ChimPipe's root directory
##############################
# It will be exported as an environmental variable since it will be used by every ChimPipe's scripts 
# to set the path to the bin, awk and bash directories. 
root=/nfs/users/rg/brodriguez/Chimeras_project/Chimeras_detection_pipeline/ChimPipe

export rootDir=$root 

# 2. Parse input arguments with getopt  
######################################
getopt=$root/bin/getopt

getoptions $0 $@ # call Function 5 and passing two parameters (name of the script and command used to call it)

# 3. Check input variables 
##########################

## Mandatory arguments
## ~~~~~~~~~~~~~~~~~~~

if [[ "$bamAsInput" != "true" ]];
then
	## A) FASTQ as input
	bamAsInput="false";
	if [[ ! -e $fastq1 ]]; then log "The mate 1 FASTQ provided does not exist. Mandatory argument --fastq_1\n" "ERROR" >&2; usagedoc; exit -1; fi
	if [[ ! -e $fastq2 ]]; then log "The mate 2 FASTQ provided does not exist. Mandatory argument --fastq_2\n" "ERROR" >&2; usagedoc; exit -1; fi
	if [[ ! -e $transcriptomeIndex ]]; then log "The transcriptome index provided does not exist. Mandatory argument -t|--transcriptome-index\n" "ERROR" >&2; usagedoc; exit -1; fi
	if [[ ! -e $transcriptomeKeys ]]; then log "The transcriptome keys provided do not exist. Mandatory argument -k|--transcriptome-keys\n" "ERROR" >&2; usagedoc; exit -1; fi
else
	## B) BAM as input
	if [[ ! -e $bam ]]; then log "The BAM provided do not exist. Mandatory argument --bam\n" "ERROR" >&2; usagedoc; exit -1; fi
fi
	
## Common
if [[ ! -e $genomeIndex ]]; then log "The genome index provided does not exist. Mandatory argument -g|--genome-index\n" "ERROR" >&2; usagedoc; exit -1; fi
if [[ ! -e $annot ]]; then log "The annotation provided does not exist. Mandatory argument -a|--annotation\n" "ERROR" >&2; usagedoc; exit -1; fi
if [[ "$lid" == "" ]]; then log "Please provide a sample identifier. Mandatory argument --sample-id\n" "ERROR" >&2; usagedoc; exit -1; fi


## Optional arguments
## ~~~~~~~~~~~~~~~~~~

# General
# ~~~~~~~

# Log level
if [[ "$logLevel" == "" ]]; 
then 
	logLevel='warn'; 
else	
	if [[ "$logLevel" != @(error|warn|info|debug) ]];
	then
		log "Please specify a proper log status [error||warn||info||debug]. Option -l|--log\n" "ERROR" >&2;
		usagedoc;
		exit -1; 
	fi
fi

# Number of threads
if [[ "$threads" == "" ]]; 
then 
	threads='1'; 
else
	if [[ ! "$threads" =~ ^[0-9]+$ ]]; 
	then
		log "Please specify a proper threading value. Option -t|--threads\n" "ERROR" >&2;
		usagedoc;
		exit -1; 
	fi
fi

if [[ "$threads" == "1" ]]; 
then 
	hthreads='1'; 
else 
	hthreads=$((threads/2)); 
fi	

# Output directory
if [[ "$outDir" == "" ]]; 
then 
	outDir=${SGE_O_WORKDIR-$PWD};
else
	if [[ ! -e "$outDir" ]]; 
	then
		log "Your output directory does not exist. Option -o|--output-dir\n" "ERROR" >&2;
		usagedoc; 
		exit -1; 
	fi	
fi

# Temporary directory
if [[ "$TMPDIR" == "" ]]; 
then 
	TMPDIR='/tmp'; 
else	
	if [[ ! -e "$TMPDIR" ]]; 
	then
		log "Your temporary directory does not exist. Option --tmp-dir\n" "ERROR" >&2;
		usagedoc; 
		exit -1; 
	fi
fi

# Clean up
if [[ "$cleanup" != "false" ]]; 
then 
    cleanup='true'; 
fi	

# Reads information:
# ~~~~~~~~~~~~~~~~~~
# Maximum read length

if [[ "$maxReadLength" == "" ]]; 
then 
    maxReadLength='150'; 
else
    if [[ ! "$maxReadLength" =~ ^[0-9]+$ ]]; 
    then
		log "Please specify a proper maximum read length value for mapping. Option --max-read-length\n" "ERROR" >&2;
		usagelongdoc;
	exit -1; 
    fi
fi

# Sequencing library type
if [[ "$readDirectionality" != "" ]];
then
    if [[ "$readDirectionality" == @(MATE1_SENSE|MATE2_SENSE) ]];
    then
		stranded=1;
    elif [[ "$readDirectionality" == "UNSTRANDED" ]];
    then
		stranded=0;
    else
		log "Please specify a proper sequencing library [UNSTRANDED|MATE1_SENSE|MATE2_SENSE]\n" "ERROR" >&2;
		usagedoc; 
	exit -1;	
    fi
else
    readDirectionality="NOT DEFINED"
fi 	

# First mapping parameters:
# ~~~~~~~~~~~~~~~~~~~~~~~~~

# Consensus splice sites for the segmental mapping
if [[ "$spliceSitesFM" == "" ]]; 
then 
    spliceSitesFM="GT+AG,GC+AG,ATATC+A.,GTATC+AT"; 
else			
    if [[ ! "$spliceSitesFM" =~ ^([ACGT.]+\+[ACGT.]+,)*([ACGT.]+\+[ACGT.]+)$ ]];
    then
	log "Please specify a proper consensus splice site sequence for the first segmental mapping. Option -C|--consensus-ss-fm\n" "ERROR" >&2;
	usagelongdoc;
	exit -1; 
    fi
fi

# Minimum split size for the segmental mapping
if [[ "$splitSizeFM" == "" ]];
then 
    splitSizeFM='15'; 
else
    if [[ ! "$splitSizeFM" =~ ^[0-9]+$ ]]; 
    then
	log "Please specify a proper minimum split size for the first segmental mapping step. Option -S|--min-split-size-fm\n" "ERROR" >&2;
	usagelongdoc; 
	exit -1; 
    fi
fi

# Refinement size for the segmental mapping
if [[ "$refinementFM" == "" ]];
then 
    refinementFM='2'; 
else
    if [[ ! "$refinementFM" =~ ^[0-9]+$ ]]; 
    then
	log "Please specify a proper refinement size for the first segmental mapping step. Option --refinement-step-size-fm\n" "ERROR" >&2;
	usagelongdoc; 
	exit -1; 
    fi
fi

# First mapping statistics
if [[ "$mapStats" == "false" ]]; 
then 
    stats="--no-stats"; 
else
	mapStats="true";
fi

# Second mapping parameters:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~

# Remap the unmapped reads from the first mappping step
if [[ "$remapUnmapped" == "false" ]]; 
then 
	extractUnmappedSM="0";	
else			
	remapUnmapped="true"; 
	extractUnmappedSM="1";
fi

# Remap multimapped reads with more than X matches from the first mappping step
if [[ "$remapMultimapped" == "true" ]]; 
then 
	extractMultimappedSM="1";
	if [[ "$nbMultimaps2Remap" == "ALL" ]];
	then
		nbMultimaps2Remap="";	
	fi
else
	remapMultimapped="false";
	extractMultimappedSM="0";
fi
	    	
# Remap uniquely mapped reads with more than X mismatches from the first mappping step
if [[ "$remapUnique" == "true" ]]; 
then 
	extractUniqueFM="1";
	extractUniqueSM="1";
else
	remapUnique="false";
	extractUniqueFM="1";
	extractUniqueSM="0";
fi

# Consensus splice sites for the segmental mapping
if [[ "$spliceSitesSM" == "" ]]; 
then 
    spliceSitesSM="GT+AG"; 
else			
    if [[ ! "$spliceSitesSM" =~ ^([ACGT.]+\+[ACGT.]+,)*([ACGT.]+\+[ACGT.]+)$ ]];
    then
	log "Please specify a proper consensus splice site sequence for the second segmental mapping. Option -c|--consensus-ss-sm\n" "ERROR" >&2;
	usagelongdoc; 
	exit -1; 
    fi
fi

# Minimum split size for the segmental mapping
if [[ "$splitSizeSM" == "" ]];
then 
    splitSizeSM='15'; 
else
    if [[ ! "$splitSizeSM" =~ ^[0-9]+$ ]]; 
    then
	log "Please specify a proper minimum split size for the second segmental mapping step. Option -s|--min-split-size-sm\n" "ERROR" >&2;
	usagelongdoc; 
	exit -1; 
    fi
fi

# Refinement size for the segmental mapping
if [[ "$refinementSM" == "" ]];
then 
    refinementSM='2'; 
else
    if [[ ! "$refinementSM" =~ ^[0-9]+$ ]]; 
    then
	log "Please specify a proper refinement size for the second segmental mapping step. Option --refinement-step-size-sm\n" "ERROR" >&2;
	usagelongdoc; 
	exit -1; 
    fi
fi

# Chimera detection phase parameters:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Consider multimapped reads for chimera detection:

if [[ $multimapped2ChimDetection == "true" ]]
then
	extractMultimappedFM="1";
	
	if [[ "$nbMultimaps2ChimDetect" == "ALL" ]];
	then
		nbMultimaps2ChimDetect="";	
	fi
else
	multimapped2ChimDetection="false";
	extractMultimappedFM="0";
fi	
	
# Filtering module configuration	
if [[ "$filterConf" == "" ]]; 
then 			
    filterConf="5,0,80,30;1,1,80,30;";		# Default
else
    if [[ ! "$filterConf" =~ ^([0-9]+,[0-9]+,[0-9]{,3},[0-9]+;){1,2}$ ]]; 
    then
	log "Please check your filtering module configuration. Option --filter-chimeras\n" "ERROR" >&2; 
	usagelongdoc; 
	exit -1;
    fi
fi

# Similarity between gene pairs file
if [[ "$simGnPairs" == "" ]]
then
    simGnPairs="NOT PROVIDED";
elif [ ! -s "$simGnPairs" ]
then 
    log "Your text file containing similarity information between gene pairs in the annotation does not exist. Option --similarity-gene-pairs\n" "ERROR" >&2; 
    usagelongdoc; 
    exit -1; 
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
#####################
# Bash 
qual=$bashDir/detect.fq.qual.sh
addXS=$bashDir/sam2cufflinks.sh
infer_library=$bashDir/infer_library_type.sh
chim1=$bashDir/find_exon_exon_connections_from_splitmappings.sh
chim2=$bashDir/find_chimeric_junctions_from_exon_to_exon_connections.sh
findGeneConnections=$bashDir/find_gene_to_gene_connections_from_pe_rnaseq.sh
sim=$bashDir/similarity_bt_gnpairs.sh

# Bin 
gemtools=$binDir/gemtools-1.7.1-i3/gemtools
gemrnatools=$binDir/gemtools-1.7.1-i3/gem-rna-tools
gtfilter=$binDir/gemtools-1.7.1-i3/gt.filter
gem2sam=$binDir/gemtools-1.7.1-i3/gem-2-sam
gt_filter=$binDir/gemtools-1.7.1-i3/gt.filter.remove
pigz=$binDir/pigz

# Awk 
correctNMfield=$awkDir/correctNMfield_SAMFromGEM.awk
SAMfilter=$awkDir/SAMfilter.awk 
addMateInfoSam=$awkDir/add_mateInfo_SAM.awk
gff2Gff=$awkDir/gff2gff.awk
bed2bedPE=$awkDir/bed2bedPE.awk
bedPECorrectStrand=$awkDir/bedPECorrectStrand.awk
bedPE2gff=$awkDir/bedPE2gff.awk
gemCorrectStrand=$awkDir/gemCorrectStrand.awk
gemToGff=$awkDir/gemsplit2gff_unique.awk
addPEinfo=$awkDir/add_PE_info.awk
AddSimGnPairs=$awkDir/add_sim_bt_gnPairs.awk
juncFilter=$awkDir/chimjunc_filter.awk


## DISPLAY PIPELINE CONFIGURATION  
##################################
printf "\n"
header="PIPELINE CONFIGURATION FOR $lid"
echo $header
eval "for i in {1..${#header}};do printf \"-\";done"
printf "\n\n"
printf "  %-34s %s\n\n" "ChimPipe Version $version"
printf "  %-34s %s\n" "***** MANDATORY ARGUMENTS *****"

if [[ "$bamAsInput" == "false" ]];
then
	## A) FASTQ as input
	printf "  %-34s %s\n" "fastq_1:" "$fastq1"
	printf "  %-34s %s\n" "fastq_2:" "$fastq2"
	printf "  %-34s %s\n" "genome-index:" "$genomeIndex"
	printf "  %-34s %s\n" "annotation:" "$annot"
	printf "  %-34s %s\n" "transcriptome-index:" "$transcriptomeIndex"
	printf "  %-34s %s\n" "transcriptome-keys:" "$transcriptomeKeys"
	printf "  %-34s %s\n\n" "sample-id:" "$lid"

	printf "  %-34s %s\n" "** Reads information **"
	printf "  %-34s %s\n" "seq-library:" "$readDirectionality"
	printf "  %-34s %s\n\n" "max-read-length:" "$maxReadLength"

	printf "  %-34s %s\n" "***** MAPPING PHASE *****"
	printf "  %-34s %s\n" "** 1st mapping **"
	printf "  %-34s %s\n" "consensus-ss-fm:" "$spliceSitesFM"
	printf "  %-34s %s\n" "min-split-size-fm:" "$splitSizeFM"
	printf "  %-34s %s\n" "refinement-step-size-fm (0:disabled):" "$refinementFM"
	printf "  %-34s %s\n\n" "stats:" "$mapStats"
else
	## B) BAM as input
	printf "  %-34s %s\n" "bam:" "$bam"
	printf "  %-34s %s\n" "genome-index:" "$genomeIndex"
	printf "  %-34s %s\n" "annotation:" "$annot"
	printf "  %-34s %s\n\n" "sample-id:" "$lid"

	printf "  %-34s %s\n" "** Reads information **"
	printf "  %-34s %s\n" "seq-library:" "$readDirectionality"
	printf "  %-34s %s\n\n" "max-read-length:" "$maxReadLength"
	
	printf "  %-34s %s\n" "***** MAPPING PHASE *****"
fi

printf "  %-34s %s\n" "** 2nd mapping **"
printf "  %-34s %s\n" "remap-unmapped:" "$remapUnmapped"
printf "  %-34s %s\n" "remap-multimapped:" "$remapMultimapped"
if [[ "$remapMultimapped" == "true" ]]; 
then 
	if [[ "$nbMultimaps2Remap" == "" ]]; 
	then
		printf "  %-34s %s\n" "(remap all the multimapped reads)"
	else
		printf "  %-34s %s\n" "(remap multimapped reads with more than $nbMultimaps2Remap matches)"
	fi
fi 

printf "  %-34s %s\n" "remap-unique:" "$remapUnique"
if [[ "$remapUnique" == "true" ]]; 
then 
	printf "  %-34s %s\n" "(remap unique mapping with more than $nbMism mismatches. Mismatches in bad-quality bases also considered)"
fi 
printf "  %-34s %s\n" "consensus-ss-fm:" "$spliceSitesSM"
printf "  %-34s %s\n" "min-split-size-fm:" "$splitSizeSM"
printf "  %-34s %s\n\n" "refinement-step-size-fm (0:disabled):" "$refinementSM"

printf "  %-34s %s\n" "***** CHIMERA DETECTION PHASE *****"
printf "  %-34s %s\n" "consider-multimapped:" "$multimapped2ChimDetection"
if [[ "$multimapped2ChimDetection" == "true" ]]; 
then 
	if [[ "$nbMultimaps2ChimDetect" == "" ]]; 
	then
		printf "  %-34s %s\n" "(consider all the multimapped reads for chimera detection)"
	else
		printf "  %-34s %s\n" "(Consider multimapped reads with $nbMultimaps2ChimDetect or less matches for chimera detection)"
	fi
fi 

printf "  %-34s %s\n" "filter-chimeras:" "$filterConf"
printf "  %-34s %s\n\n" "similarity-gene-pairs:" "$simGnPairs"

printf "  %-34s %s\n" "***** GENERAL *****"
printf "  %-34s %s\n" "output-dir:" "$outDir"
printf "  %-34s %s\n" "tmp-dir:" "$TMPDIR"
printf "  %-34s %s\n" "threads:" "$threads"
printf "  %-34s %s\n" "log:" "$logLevel"
printf "  %-34s %s\n\n" "no-cleanup:" "$cleanup"

	
## START CHIMPIPE
#################
header="Executing ChimPipe $version for $lid"
echo $header
eval "for i in {1..${#header}};do printf \"-\";done"
printf "\n\n"
pipelineStart=$(date +%s)

# 0) Preliminary steps
######################

if [[ "$bamAsInput" == "false" ]];
then
	step="PRELIM"
	log "Determining the offset quality of the reads for ${lid}..." $step
	run "quality=\`$qual $fastq1 | awk '{print \$2}'\`" "$ECHO" 
	log " The read quality is $quality\n"
	log "done\n"
else
	quality="33"
fi

b=`basename $annot`
b2tmp=${b%.gtf}
b2=${b2tmp%.gff}
    	
    	
# 1) First mapping step. Map all the reads to the genome, to the transcriptome and de-novo, using the 
#################################################################################################
# gemtools RNA-Seq pipeline but with max intron size larger than the biggest chromosome, and with 
##################################################################################################  
# a number of mismatches of round(read_length/6) an edit distance of round(read_length/20) 
##########################################################################################
# outputs are: 
##############
# - $outDir/${lid}.map.gz 
# - $outDir/${lid}_raw_chrSorted.bam

if [[ "$bamAsInput" == "false" ]];
then
	bamFirstMap=$outDir/${lid}_firstMap.bam
	if [ ! -s $bamFirstMap ]; 
	then
		runGemtoolsRnaPipeline		## Call function to run the gemtools rna-pipeline
	else
    	printHeader "First mapping BAM file already exists... skipping first mapping step"
	fi
else
	printHeader "BAM file provided as input... skipping first mapping step"
	bamFirstMap=$bam	
fi

### Check in which fields are the number of mappings and the number of mismatches 

# Comment: samtools view -F 4 ** filter out unmapped reads 

NHfield=`samtools view -F 4 $bamFirstMap | head -1 | awk 'BEGIN{field=1;}{while(field<=NF){ if ($field ~ "NH:i:"){print field;} field++}}'`
NMfield=`samtools view -F 4 $bamFirstMap | head -1 | awk 'BEGIN{field=1;}{while(field<=NF){ if ($field ~ "NM:i:"){print field;} field++}}'`


# 2) Infer the sequencing library protocol used (UNSTRANDED, MATE2_SENSE OR MATE1_SENSE) 
########################################################################################
# from a subset with the 1% of the mapped reads. 
#################################################
# Outputs are: 
##############
# - variables $readDirectionality and $stranded

if [[ "$readDirectionality" == "NOT DEFINED" ]]; 
then 
    step="INFER-LIBRARY"
    startTimeLibrary=$(date +%s)
    printHeader "Executing infer library type step" 
    log "Infering the sequencing library protocol from a random subset with 1 percent of the mapped reads..." $step
    read fraction1 fraction2 other <<<$(bash $infer_library $bamFirstMap $annot)
    log "done\n"
    log "Fraction of reads explained by 1++,1--,2+-,2-+: $fraction1\n" $step
    log "Fraction of reads explained by 1+-,1-+,2++,2--: $fraction2\n" $step
    log "Fraction of reads explained by other combinations: $other\n" $step 
    
    # Turn the percentages into integers
    fraction1_int=${fraction1/\.*};
    fraction2_int=${fraction2/\.*};
    other_int=${other/\.*};
    
    # Infer the sequencing library from the mapping distribution. 
    if [ "$fraction1_int" -ge 70 ]; # MATE1_SENSE protocol
    then 
		readDirectionality="MATE1_SENSE";
		stranded=1;
		echo $readDirectionality;	
    elif [ "$fraction2_int" -ge 70 ];
    then	
		readDirectionality="MATE2_SENSE"; # MATE2_SENSE protocol
		stranded=1;
    elif [ "$fraction1_int" -ge 40 ] && [ "$fraction1_int" -le 60 ];
    then
		if [ "$fraction2_int" -ge 40 ] && [ "$fraction2_int" -le 60 ]; # UNSTRANDED prototol
		then
	    	readDirectionality="UNSTRANDED";
	    	stranded=0;
		else
	    	log "ChimPipe is not able to determine the library type. Ask your data provider and use the option -l|--seq-library\n" "ERROR" >&2;
	    	usagelongdoc
	    	exit -1	
		fi
    else
		log "ChimPipe is not able to determine the library type. Ask your data provider and use the option -l|--seq-library\n" "ERROR" >&2;
		usagelongdoc
		exit -1	
    fi
    log "Sequencing library type: $readDirectionality\n" $step 
    log "Strand aware protocol (1: yes, 0: no): $stranded\n" $step 
    endTimeLibrary=$(date +%s)
    printHeader "Sequencing library inference for $lid completed in $(echo "($endTimeLibrary-$startTimeLibrary)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Sequencing library type provided by the user...skipping library inference step"
fi

# 3) Produce a filtered BAM file with the alignments of the reads which will not be remapped. 
#############################################################################################
# Add the XS optional field to specify the read directionality for cufflinks
##############################################################################
# Output is: 
############
# - $outDir/${lid}_filtered_chrSorted.bam

## Comment: samtools view -F 256 ** filter out secondary alignments (multimapped reads represented as one primary alignment plus several secondary alignments) 

filteredBam=$outDir/${lid}_filtered_firstMap.bam
		
if [ ! -s $filteredBam ];
then
	step="BAM-FILTERING"
    startTime=$(date +%s)
    printHeader "Executing bam filtering step"
	log "Produce a filtered BAM file with the mappings of the reads which will not be remapped..." $step
	run "samtools view -h -F 256 $bamFirstMap | awk -v OFS="'\\\t'" -f $correctNMfield | awk -v OFS="'\\\t'" -v higherThan="0" -v unmapped="0" -v multimapped=$extractMultimappedFM -v nbMatches=$nbMultimaps2ChimDetect -v unique=$extractUniqueFM -v nbMism=$nbMism -v NHfield=$NHfield -v NMfield=$NMfield -f $SAMfilter | $addXS $readDirectionality | samtools view -@ $threads -Sb - | samtools sort -@ $threads -m 4G - ${filteredBam%.bam}" "$ECHO"
	log "done\n"
	if [ -s $filteredBam ]; 
   	then
       	log "Computing md5sum for bam file..." $step
       	run "md5sum $filteredBam > $filteredBam.md5" "$ECHO"
       	log "done\n"
   	else
       	log "Error filtering the bam file\n" "ERROR" 
       	exit -1
   	fi
   	endTime=$(date +%s)
   	printHeader "Bam filtering step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
   	printHeader "Filtered bam file already exists...skipping conversion step"
fi


# 4) Extract reads from the raw BAM produced in the first mapping for a second 
###############################################################################
# split-mapping attemp allowing split-mappings in different chromosomes, strands 
#################################################################################
# and non genomic order. Produce a FASTQ file with them. 
########################################################
# output is: 
############
# - $outDir/${lid}_reads2remap.fastq

## Comment: samtools view -F 256 ** filter out secondary alignments (multimapped reads represented as one primary alignment plus several secondary alignments)

reads2remap=$outDir/${lid}_reads2remap.fastq
		
if [ ! -s $reads2remap ]; 
then
	step="READS2REMAP"
	startTime=$(date +%s)
	printHeader "Executing extract reads to remap step" 
	log "Extracting reads from the raw BAM for a second split-mapping step..." $step
	run "samtools view -h -F 256 $bamFirstMap | awk -v OFS="'\\\t'" -f $correctNMfield | awk -v OFS="'\\\t'" -v higherThan="1" -v unmapped=$extractUnmappedSM -v multimapped=$extractMultimappedSM -v nbMatches=$nbMultimaps2Remap -v unique=$extractUniqueSM -v nbMism=$nbMism -v NHfield=$NHfield -v NMfield=$NMfield -f $SAMfilter | awk -v OFS="'\\\t'" -f $addMateInfoSam | samtools view -@ $threads -bS - | bedtools bamtofastq -i - -fq $reads2remap" "$ECHO"
	log "done\n" 
    if [ -s $reads2remap ]; 
    then
    	log "Computing md5sum for the fastq file..." $step
    	run "md5sum $reads2remap > $reads2remap.md5" "$ECHO"
        log "done\n"
    else
        log "Error extracting the reads\n" "ERROR" 
        exit -1
    fi
    endTime=$(date +%s)
	printHeader "Extracting reads completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "FASTQ file with reads to remap already exists... skipping extracting reads to remap step"
fi

# 5) Second split-mapping attemp. Remap the extracted reads allowing reads 
###########################################################################
# to split in different chromosomes, strands and non genomic order.
###################################################################
# output is: 
############
# - $outDir/SecondMapping/${lid}.remapped.map

gemSecondMap=$outDir/SecondMapping/${lid}_secondMap.map

if [ ! -s $gemSecondMap ];
then
	step="SECOND-MAP"
	startTime=$(date +%s)
	printHeader "Executing second split-mapping step"
	log "Remapping reads allowing them to split-map in different chromosomes, strand and non genomic order..." $step	
	run "$gemrnatools split-mapper -I $genomeIndex -i $reads2remap -q 'offset-$quality' -o ${gemSecondMap%.map} -t 10 -T $threads --min-split-size $splitSizeSM --refinement-step-size $refinementSM --splice-consensus $spliceSitesSM  > $outDir/SecondMapping/$lid.gem-rna-mapper.out" "$ECHO"
	log "done\n" 
	if [ -s $gemSecondMap ]; 
	then
    	log "Computing md5sum for the gem file with the second mappings..." $step
    	run "md5sum $gemSecondMap > $gemSecondMap.md5" "$ECHO"
        log "done\n"
    else
        log "Error in the second mapping\n" "ERROR" 
        exit -1
    fi
    endTime=$(date +%s)
	printHeader "Unmapped reads mapping step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
	printHeader "Second mapping GEM file already exists... skipping extracting second mapping step"
fi
	
# 6) extract the reads mapping both uniquely and in 2 blocks from the bam file of "normal" mappings and convert in gff.gz
#########################################################################################################################
# output is: 
############
# - $outDir/FromFirstBam/${lid}_splitmappings_2blocks_firstMap.gff.gz

gffFromBam=$outDir/FromFirstBam/${lid}_splitmappings_2blocks_firstMap.gff.gz

if [ ! -s $gffFromBam ]; 
then
	step="FIRST-CONVERT"
	startTime=$(date +%s)
	printHeader "Executing conversion of the bam into gff step"
	log "Generating a ".gff.gz" file from the normal mappings containing the reads split-mapping both uniquely and in 2 blocks..." $step
	bedtools bamtobed -i $filteredBam -bed12 | awk '$10==2' | awk -v rev='1' -f $bed2bedPE | awk -v readDirectionality=$readDirectionality  -f $bedPECorrectStrand | awk -f $bedPE2gff | awk -f $gff2Gff | gzip > $gffFromBam 
	log "done\n"
	if [ -s $gffFromBam ]; 
	then
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
	printHeader "Gff from first mapping BAM file already exists... skipping conversion step"
fi

# 7) extract the reads mapping both uniquely and in 2 blocks from the map file of "atypical" mappings and convert in gff.gz
#########################################################################################################################
# output is: 
############
# - $outDir/FromSecondMapping/${lid}_splitmappings_2blocks_secondMap.gff.gz

gffFromMap=$outDir/FromSecondMapping/${lid}_splitmappings_2blocks_secondMap.gff.gz

if [ ! -s $gffFromMap ]; 
then
	step="SECOND-CONVERT"
	startTime=$(date +%s)
	printHeader "Executing conversion of the gem into gff step"
	log "Generating a ".gff.gz" file from the atypical mappings containing the reads split-mapping both uniquely and in 2 blocks..." $step	
	run "awk -v readDirectionality=$readDirectionality -f $gemCorrectStrand $gemSecondMap | awk -v rev="0" -f $gemToGff | awk -f $gff2Gff | gzip > $gffFromMap" "$ECHO"
	log "done\n" 
	if [ -s $gffFromMap ]; 
	then
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
	printHeader "Gff from second mapping MAP file already exists... skipping conversion step"
fi


# 8) put the path to the "normal" and "atypical" gff.gz files in a same txt file for chimsplice 
#############################################################################################
paths2chimsplice=$outDir/split_mapping_file_sample_$lid.txt

run "echo $gffFromBam > $paths2chimsplice" "$ECHO"
run "echo $gffFromMap >> $paths2chimsplice" "$ECHO"

# 9) run chimsplice on 4) and 5)
###############################
# - $outDir/Chimsplice/chimeric_junctions_report_$lid.txt
# - $outDir/Chimsplice/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt
# - $outDir/Chimsplice/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn_morethan10staggered.txt

exonConnections1=$outDir/Chimsplice/exonA_exonB_with_splitmapping_part1overA_part2overB_readlist_sm1list_sm2list_staggeredlist_totalist_${lid}_splitmappings_2blocks_firstMap.txt.gz
exonConnections2=$outDir/Chimsplice/exonA_exonB_with_splitmapping_part1overA_part2overB_readlist_sm1list_sm2list_staggeredlist_totalist_${lid}_splitmappings_2blocks_secondMap.txt.gz

printHeader "Executing Chimsplice step"
if [ ! -s $exonConnections1 ] || [ ! -s $exonConnections2 ]; 
then
	step="CHIMSPLICE"
	startTime=$(date +%s)
	log "Finding exon to exon connections from the ".gff.gz" files containing the "normal" and "atypical" mappings..." $step
	run "$chim1 $paths2chimsplice $annot $outDir/Chimsplice $stranded" "$ECHO"
	log "done\n" 
	if [ ! -s $exonConnections1 ] || [ ! -s $exonConnections2 ]; 
	then
        log "Error running chimsplice\n" "ERROR" 
        exit -1
    fi
    endTime=$(date +%s)
	printHeader "Find exon to exon connections step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
	printHeader "Exon to exon connections file already exists... skipping step"
fi

chimJunctions=$outDir/Chimsplice/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt
if [ ! -s $chimJunctions ]; 
then
	step="CHIMSPLICE"
	startTime=$(date +%s)
	log "Finding chimeric junctions from exon to exon connections..." $step
	run "$chim2 $paths2chimsplice $genomeIndex $annot $outDir/Chimsplice $stranded $spliceSitesFM > $outDir/Chimsplice/chimeric_junctions_report_$lid.txt 2> $outDir/Chimsplice/find_chimeric_junctions_from_exon_to_exon_connections_$lid.err" "$ECHO"
	log "done\n" 
	if [ ! -s $chimJunctions ]; 
	then
        log "Error running chimsplice\n" "ERROR" 
        exit -1
    fi
    endTime=$(date +%s)
	printHeader "Find chimeric junctions step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
	printHeader "Chimeric Junctions file already exists... skipping step"
fi

# 9.1) Find gene to gene connections supported by paired-end mappings from the bam file of "normal" mappings with the number of mappings supporting the connection.  
#################################################################################################################################################################
# For a connection g1 to g2 to exist there must be at least one mapping where the first mate is strandedly (if data is stranded) overlapping with an exon 
#########################################################################################################################################################
# of g1 and the second mate is (strandedly if data is stranded) overlapping with an exon of g2
##############################################################################################
# - $outDir/PE/readid_gnlist_whoseexoverread_noredund.txt.gz
# - $outDir/PE/readid_twomateswithgnlist_alldiffgnpairs_where_1stassociatedto1stmate_and2ndto2ndmate.txt.gz
# - $outDir/PE/pairs_of_diff_gn_supported_by_pereads_nbpereads.txt

PEsupport=$outDir/PE/pairs_of_diff_gn_supported_by_pereads_nbpereads.txt

printHeader "Executing find gene to gene connections from PE mappings step"
if [ ! -s $PEsupport ];
then
	step="PAIRED-END"
	startTime=$(date +%s)
	log "Finding gene to gene connections supported by paired-end mappings from the ".bam" containing reads mapping in a unique and continuous way..." $step
	run "$findGeneConnections $filteredBam $annot $outDir/PE $readDirectionality" "$ECHO"
	log "done\n" 
	if [ ! -s $PEsupport ]; 
	then
        log "Error finding gene to gene connections\n" "ERROR" 
        exit -1
    fi
    endTime=$(date +%s)
	printHeader "Find gene to gene connections from PE mappings step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
	printHeader "Gene to gene connections file already exists... skipping step"
fi

# 9.2) Add gene to gene connections information to chimeric junctions matrix
##########################################################################
# - $outDir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_PEinfo_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt

chimJunctionsPE=$outDir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_PEinfo_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt

if [ ! -s $chimJunctionsPE ];
then
	step="PAIRED-END"
	startTime=$(date +%s)
	log "Adding PE information to the matrix containing chimeric junction candidates..." $step
	run "awk -v fileRef=$PEsupport -f $addPEinfo $chimJunctions 1> $outDir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_PEinfo_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt
" "$ECHO"
	log "done\n" 
	if [ ! -s $chimJunctionsPE ]; then
        log "Error adding PE information\n" "ERROR" 
        exit -1
    fi
    endTime=$(date +%s)
	printHeader "Add PE information step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
	printHeader "Chimeric junction matrix with PE information already exists... skipping step"
fi

# 10) Compute the gene similarity matrix in case the user does not provide it
############################################################################
if [ ! -e "$simGnPairs" ]
then
    step="PRE-SIM"
    startTime=$(date +%s)
    log "Computing similarity between annotated genes..." $step
    run "$sim $annot $genomeIndex 4" "$ECHO"
    log "done\n" 			
    endTime=$(date +%s)
    simGnPairs=$outDir/$b2\_gene1_gene2_alphaorder_pcentsim_lgalign_trpair.txt
    printHeader "Computing similarity between annotated gene step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
fi

# 11) Add information regarding the sequence similarity between connected genes
##############################################################################
# - $outDir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_PEinfo_maxLgalSim_maxLgal_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt

chimJunctionsSim=$outDir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_PEinfo_maxLgalSim_maxLgal_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt

if [  ! -e "$chimJunctionsSim" ]
then		
    if [ -e "$simGnPairs" ]
    then
	step="SIM"
	startTime=$(date +%s)
	log "Adding sequence similarity between connected genes information to the chimeric junction matrix..." $step
	run "awk -v fileRef=$simGnPairs -f $AddSimGnPairs $chimJunctionsPE > $chimJunctionsSim" "$ECHO"
	log "done\n" 
	if [ ! -s $chimJunctionsSim ]
	then
	    log "Error adding similarity information, file is empty\n" "ERROR" 
	    exit -1
	fi
	endTime=$(date +%s)
	printHeader "Add sequence similarity information step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
    else 
	printHeader "Similarity information between the gene pairs in the annotation is not provided... skipping step"
    fi
else
    printHeader "Chimeric junction matrix with similarity information already exists... skipping step"
fi


# 12) Produce a matrix containing chimeric junction candidates with a header in the first row
#############################################################################################
# - $outDir/chimeric_junctions_candidates.txt

chimJunctionsCandidates=$outDir/chimeric_junctions_candidates_${lid}.txt

if [ ! -s "$chimJunctionsCandidates" ]
then		
    step="HEADER"
    log "Adding a header to the matrix containing the chimeric junction candidates..." $step
    if [ -s "$chimJunctionsSim" ]
    then		
		run "awk 'BEGIN{print \"juncId\", \"nbstag\", \"nbtotal\", \"maxbeg\", \"maxEnd\", \"samechr\", \"samestr\", \"dist\", \"ss1\", \"ss2\", \"gnlist1\", \"gnlist2\", \"gnname1\", \"gnname2\", \"bt1\", \"bt2\", \"PEsupport\", \"maxSim\", \"maxLgal\";}{print \$0;}' $chimJunctionsSim 1> $chimJunctionsCandidates" "$ECHO"
	log "done\n"
    else
		if [ -s "$chimJunctionsPE" ]
		then
	   		run "awk 'BEGIN{print \"juncId\", \"nbstag\", \"nbtotal\", \"maxbeg\", \"maxEnd\", \"samechr\", \"samestr\", \"dist\", \"ss1\", \"ss2\", \"gnlist1\", \"gnlist2\", \"gnname1\", \"gnname2\", \"bt1\", \"bt2\", \"PEsupport\";}{print \$0;}' $chimJunctionsPE 1> $chimJunctionsCandidates" "$ECHO"
	    	log "done\n" 	
		else
	    	log "Error, intermediate file: $chimJunctionsSim or $chimJunctionsPE is missing\n" "ERROR" 
	    exit -1			
		fi
    fi 
else
    printHeader "Header already already added... skipping step"
fi

# 13) Filter out chimera candidates to produce a final set of chimeric junctions
################################################################################
# - $outDir/chimeric_junctions.txt

chimJunctions=$outDir/chimeric_junctions_${lid}.txt

if [ ! -s "$chimJunctions" ]; 
then
	step="FILTERING MODULE"
	startTime=$(date +%s)
	log "Filtering out chimera candidates to produce a final set of chimeric junctions..." $step
	awk -v filterConf=$filterConf -f $juncFilter $chimJunctionsCandidates > $chimJunctions
	log "done\n" 
	if [ ! -s $chimJunction ]; 
	then
		log "Error filtering chimeric junction candidates\n" "ERROR" 
    	exit -1
	fi
	
	endTime=$(date +%s)
	printHeader "Filtering module step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else 
	printHeader "Chimeric junction candidates already filtered... skipping step"
fi


# 14) Clean up
###############

if [[ "$cleanup" == "true" ]]; 
then 
	step="CLEAN-UP"
	startTime=$(date +%s)
	log "Removing intermediate files..." $step
	rm -r $outDir/FromFirstBam $outDir/FromSecondMapping $outDir/Chimsplice $outDir/PE
	rm $outDir/$b2\_tr.fasta*
	rm $outDir/${lid}.junctions $outDir/${lid}_denovo-index.junctions $outDir/${lid}_denovo-index.junctions.keys $outDir/${lid}_gtf.junctions $gemFirstMapping $gemFirstMapping.md5 $filteredGem $filteredGem.md5 $unmappedReads $unmappedReads.md5 $paths2chimsplice $chimJunctionsSim $chimJunctionsCandidates $chimJunctionsPE
	log "done\n" 
	endTime=$(date +%s)
	printHeader "Clean up step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else 
	printHeader "No clean up mode... skipping step"
fi	


# 15) END
#########
pipelineEnd=$(date +%s)
printHeader "Chimera Mapping pipeline for $lid completed in $(echo "($pipelineEnd-$pipelineStart)/60" | bc -l | xargs printf "%.2f\n") min "

# disable extglob
shopt -u extglob

exit 0

