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

*** MANDATORY 
		
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
	--sample-id			<STRING>	Sample identifier (the output files are named according to this id).  
	
*** [OPTIONS] can be:

* General: 
	--threads			<INTEGER>	Number of threads to use. Default 1.
	-o|--output-dir			<PATH>		Output directory. Default current working directory. 
	--tmp-dir			<PATH>		Temporary directory. Default /tmp.	
	--no-cleanup					Keep intermediate files. 		
	-h|--help					Display partial usage information, only mandatory plus general arguments.
	-f|--full-help					Display full usage information with additional options. 

help
}


# Function 2. Print all the other options
##########################################
function usage_long
{
cat <<help
* Reads information:
	--max-read-length		<INTEGER>	Maximum read length. This is used to create the de-novo transcriptome and acts as an upper bound. Default 150.
	-l|--seq-library 		<STRING> 	Type of sequencing library [MATE1_SENSE | MATE2_SENSE | UNSTRANDED].
                        				UNSTRANDED for not strand-specific protocol (unstranded data) and the others for the different types 
							of strand-specific protocols (stranded data).
* Mapping phase 
                        				
  First mapping:
	-C|--consensus-ss-fm		<(COUPLE_1)>, ... ,<(COUPLE_s)>	with <couple> := <donor_consensus>+<acceptor_consensus>
                                 			List of couples of donor/acceptor splice site consensus sequences. Default='GT+AG,GC+AG,ATATC+A.,GTATC+AT'
	-S|--min-split-size-fm		<INTEGER>	Minimum split size for the segmental mapping steps. Default 15.
	--refinement-step-size-fm   	<INTEGER>   	If not mappings are found a second attempt is made by eroding "N" bases toward the ends of the read. 
							A value of 0 disables it. Default 2. 
	--no-stats					Disable mapping statistics. Default enabled.

  Second Mapping:
	-c|--consensus-ss-sm		<(COUPLE_1)>, ... ,<(COUPLE_s)>	List of couples of donor/acceptor splice site consensus sequences. Default='GT+AG'
	-s|--min-split-size-sm		<INTEGER>	Minimum split size for the segmental mapping steps. Default 15.
	--refinement-step-size-sm   	<INTEGER>   	If not mappings are found a second attempt is made by eroding "N" bases toward the ends of the read. 
							A value of 0 disables it. Default 2. 
    
* Chimera detection phase 

  Classification:
    --readthrough-max-dist 		<INTEGER> 	Maximum distance between donor and acceptor sites to classify a chimeric junction as readthrought. 
    							Default 100000.
	
  Filters:
	--total-support 		<INTEGER> 	Minimum number of total supporting evidences (spanning reads + consistent paired-ends). Default 3.
	--spanning-reads		<INTEGER>  	Minimum number of junction spanning reads. Default 1.
	--consistent-pairs		<INTEGER>	Minimum number of consistent paired-ends. Default 1.
	
	--total-support-novel-ss 	<INTEGER> 	Minimum number of total supporting evidences if novel splice-sites. Default 6.
	--spanning-reads-novel-ss	<INTEGER>  	Minimum number of junction spanning reads if novel splice-sites. Default 3.
	--consistent-pairs-novel-ss	<INTEGER>	Minimum number of consistent paired-ends if novel splice-sites. Default 3.
	
	--perc-staggered		<PERCENTAGE>	Minimum percentage of staggered reads. Default 0 (not enabled).
	--perc-multimappings		<PERCENTAGE>	Maximum percentage of multimapped spanning reads. Default 100 (not enabled).
	--perc-inconsistent-pairs	<PERCENTAGE>	Maximum percentage of inconsistent paired ends. Default 100 (not enabled).
	
	--similarity			<COUPLE>	Couple "maximum alignment length" / "maximum similarity percentage" for gene pair similarity based filtering. 
							Default='30+90'.   	

	--biotype			<STRING>,...,<STRING>	Black list of gene biotypes. Chimeric junctions involving genes with their biotype in the list will be discared. Default: 'pseudogene,polymorphic_pseudogene,IG_C_pseudogene,IG_J_pseudogene,IG_V_pseudogene,TR_J_pseudogene,TR_V_pseudogene'  
	
	




	
  Files:	
	--similarity-gene-pairs		<TEXT>		Text file with similarity information between the gene pairs in the annotation.
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
# - ${lid}_firstMap.stats.txt
# - ${lid}_firstMap.stats.json


function firstMapping_FASTQinput {

	# 1.1) Produce a filtered and sorted bam file with the aligments
	####################################################	
    if [ ! -s $gemFirstMap ]; 
    then
	step="FIRST-MAP"
	startTimeFirstMap=$(date +%s)
	printHeader "Executing first mapping step"    
	
		## Copy needed files to TMPDIR
	copyToTmp "index,annotation,t-index,keys"
	
	log "Running GEMtools rna pipeline on ${lid}..." $step
	run "$gemtools --loglevel $logLevel rna-pipeline -f $fastq1 $fastq2 -i $TMPDIR/`basename $genomeIndex` -a $TMPDIR/`basename $annot` -r $TMPDIR/`basename $transcriptomeIndex` -k $TMPDIR/`basename $transcriptomeKeys` -q $quality --max-read-length $maxReadLength --max-intron-length $longestChrLength --min-split-size $splitSizeFM --refinement-step $refinementFM --junction-consensus $spliceSitesFM --no-filtered --no-bam --no-xs $stats --no-count -n `basename ${gemFirstMap%.map.gz}` --compress-all --output-dir $TMPDIR -t $threads >> $firstMappingDir/${lid}_firstMap.log 2>&1" "$ECHO" 
   
	if [ -s $TMPDIR/`basename $gemFirstMap` ]; 
	then 
			# Copy files from temporary to output directory
	    run "cp $TMPDIR/`basename $gemFirstMap` $gemFirstMap" "$ECHO"	       	    
	    endTimeFirstMap=$(date +%s)        		
	    printHeader "First mapping for $lid completed in $(echo "($endTimeFirstMap-$startTimeFirstMap)/60" | bc -l | xargs printf "%.2f\n") min"
	else
       	    log "Error running the GEMtools pipeline\n" "ERROR"
       	    exit -1
	fi		
	else
    	printHeader "First mapping GEM file already exists... skipping first mapping step"
	fi
		
	# 1.2) Produce a filtered GEM file with the aligments
	####################################################
	gemFirstMapFiltered=$firstMappingDir/${lid}_firstMap_filtered.map.gz
		
	if [ ! -s $gemFirstMapFiltered ];
	then
		step="FILTER"
		startTime=$(date +%s)
		printSubHeader "Executing filtering GEM alignments step"
		
		## Compute the maximum number of mismatches for mapping filtering 
		# based on the average read-pair length and a given percentage of mismatches
		read nbMism <<<$(bash $nbMismLen $fastq1 $fastq2 '4')
		
		## Apply the filtering
		run "$gtFilterRemove -i $gemFirstMap --max-matches 10 --max-levenshtein-error $nbMism -t $hthreads | $pigz -p $hthreads > $gemFirstMapFiltered" "$ECHO"
		
		if [ -s $gemFirstMapFiltered ];
		then
			endTime=$(date +%s)
			printSubHeader "Filtering step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
		else
			log "Error doing the filtering\n" "ERROR"
			exit -1
		fi
	else
		printSubHeader "Filtered GEM file already exist... skipping filtering step"
	fi
	
	# 1.3) Convert the SAM into a sortered BAM file
	#################################################	
	step="GEM2BAM"
	startTime=$(date +%s)
	printSubHeader "Executing conversion GEM to BAM step"

	## Copy needed files to TMPDIR
	copyToTmp "index"	
	log "Converting $lid to bam..." $step
	run "$pigz -p $hthreads -dc $gemFirstMapFiltered | $gem2sam -T $hthreads -I $TMPDIR/`basename $genomeIndex` --expect-paired-end-reads -q offset-$quality -l | samtools view -@ $threads -bS - | samtools sort -@ $threads -m 4G - ${bamFirstMap%.bam} >> $firstMappingDir/${lid}_map2bam_conversion.log 2>&1" "$ECHO"
	if [ -s $bamFirstMap ];
	then
		endTime=$(date +%s)
		printSubHeader "Conversion step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
	else
		log "Error producing the bam file\n" "ERROR"
		exit -1
	fi
		
}


# Function 9. Parse user's input
################################
function getoptions {
ARGS=`$getopt -o "g:a:t:k:o:hfl:C:S:c:s:" -l "fastq_1:,fastq_2:,bam:,genome-index:,annotation:,transcriptome-index:,transcriptome-keys:,sample-id:,log:,threads:,output-dir:,tmp-dir:,no-cleanup,dry,help,full-help,max-read-length:,seq-library:,consensus-ss-fm:,min-split-size-fm:,refinement-step-size-fm:,no-stats,consensus-ss-sm:,min-split-size-sm:,refinement-step-size-sm:,readthrough-max-dist:,total-support:,spanning-reads:,consistent-pairs:,total-support-novel-ss:,spanning-reads-novel-ss:,consistent-pairs-novel-ss:,perc-staggered:,perc-multimappings:,perc-inconsistent-pairs:,similarity:,biotype:,similarity-gene-pairs:" \
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
              bamFirstMap=$2
              bamAsInput="TRUE"
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
	  cleanup="FALSE";
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
	      mapStats="FALSE"  
	  fi
	  shift;;
      
	# Second mapping parameters:
      
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
	
	# Chimera detection phase parameters 
  	# Classification:
  	  --readthrough-max-dist)
  	  if [ -n "$2" ];
	  then
	      readthroughMaxDist=$2
	  fi
	  shift 2;;
  	
  	# Filters:
  	  --total-support)
	  if [ -n "$2" ];
	  then
	      minNbTotal=$2
	  fi
	  shift 2;;
	  
	  --spanning-reads)
	  if [ -n "$2" ];
	  then
	      minNbSpanning=$2
	  fi
	  shift 2;;
	  
	  --consistent-pairs)
	  if [ -n "$2" ];
	  then
	      minNbConsistentPE=$2
	  fi
	  shift 2;;
	  
	  --total-support-novel-ss)
	  if [ -n "$2" ];
	  then
	      minNbTotalNovelSS=$2
	  fi
	  shift 2;;
	  
	  --spanning-reads-novel-ss)
	  if [ -n "$2" ];
	  then
	      minNbSpanningNovelSS=$2
	  fi
	  shift 2;;
	  
	  --consistent-pairs-novel-ss)
	  if [ -n "$2" ];
	  then
	      minNbConsistentPENovelSS=$2
	  fi
	  shift 2;;
	  
	  --perc-staggered)
	  if [ -n "$2" ];
	  then
	      minPercStaggered=$2
	  fi
	  shift 2;;
	  
	  --perc-multimappings)
	  if [ -n "$2" ];
	  then
	      maxPercMultimaps=$2
	  fi
	  shift 2;;
	  
	  --perc-inconsistent-pairs)
	  if [ -n "$2" ];
	  then
	      maxPercInconsistentPE=$2
	  fi
	  shift 2;;
	  
	  --similarity)
	  if [ -n "$2" ];
	  then
	      similarityConf=$2
	  fi
	  shift 2;;
	  
	  --biotype)
	  if [ -n "$2" ];
	  then
	      biotype=$2
	  fi
	  shift 2;;
	  	
	  	
	  # Files:	
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
version=v0.9.5

# Enable extended pattern matching 
shopt -s extglob

# 1. ChimPipe's root directory
##############################
# to set the path to the bin, awk and bash directories. 

path="`dirname \"$0\"`"              # relative path
rootDir="`( cd \"$path\" && pwd )`"  # absolute path

if [ -z "$rootDir" ] ; 
then
  # error; for some reason, the path is not accessible
  # to the script
  log "Path not accessible to the script\n" "ERROR" 
  exit 1  # fail
fi


# 2. Parse input arguments with getopt  
######################################
getopt=$rootDir/bin/getopt

getoptions $0 $@ # call Function 5 and passing two parameters (name of the script and command used to call it)

# 3. Check input variables 
##########################

## Mandatory arguments
## ~~~~~~~~~~~~~~~~~~~

if [[ "$bamAsInput" != "TRUE" ]];
then
	## A) FASTQ as input
	bamAsInput="FALSE";
	if [[ ! -e $fastq1 ]]; then log "The mate 1 FASTQ provided does not exist. Mandatory argument --fastq_1\n" "ERROR" >&2; usagedoc; exit -1; fi
	if [[ ! -e $fastq2 ]]; then log "The mate 2 FASTQ provided does not exist. Mandatory argument --fastq_2\n" "ERROR" >&2; usagedoc; exit -1; fi
	if [[ ! -e $transcriptomeIndex ]]; then log "The transcriptome index provided does not exist. Mandatory argument -t|--transcriptome-index\n" "ERROR" >&2; usagedoc; exit -1; fi
	if [[ ! -e $transcriptomeKeys ]]; then log "The transcriptome keys provided do not exist. Mandatory argument -k|--transcriptome-keys\n" "ERROR" >&2; usagedoc; exit -1; fi
else
	## B) BAM as input
	if [[ ! -e $bamFirstMap ]]; then log "The BAM provided do not exist. Mandatory argument --bam\n" "ERROR" >&2; usagedoc; exit -1; fi
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
if [[ "$cleanup" != "FALSE" ]]; 
then 
    cleanup='TRUE'; 
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
    readDirectionality="UNKNOWN"
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
if [[ "$mapStats" == "FALSE" ]]; 
then 
    stats="--no-stats"; 
else
	mapStats="TRUE";
fi

# Second mapping parameters:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~

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

### Classification:
if [[ "$readthroughMaxDist" == "" ]];
then 
    readthroughMaxDist=100000;
else
	if [[ ! "$readthroughMaxDist" =~ ^[0-9]+$ ]]; 
    then
	log "Please specify a proper readthrough maximum distance. Option --readthrough-max-dist\n" "ERROR" >&2;
	usagelongdoc; 
	exit -1; 
    fi
fi
 
### Filters configuration	
	      
# Minimum number of total supporting evidences (spanning reads + discordant PE)
if [[ "$minNbTotal" == "" ]];
then 
    minNbTotal=3;
else
	if [[ ! "$minNbTotal" =~ ^[0-9]+$ ]]; 
    then
	log "Please specify a proper minimum supporting evidence. Option --total-support\n" "ERROR" >&2;
	usagelongdoc; 
	exit -1; 
    fi
fi

# Minimum number of spanning reads
if [[ "$minNbSpanning" == "" ]];
then 
    minNbSpanning=1;
else
	if [[ ! "$minNbSpanning" =~ ^[0-9]+$ ]]; 
    then
	log "Please specify a proper minimum number of spanning reads. Option --spanning-reads\n" "ERROR" >&2;
	usagelongdoc; 
	exit -1; 
    fi
fi
	
# Minimum number of consistent paired-ends
if [[ "$minNbConsistentPE" == "" ]];
then 
    minNbConsistentPE=1;
else
	if [[ ! "$minNbConsistentPE" =~ ^[0-9]+$ ]]; 
    then
	log "Please specify a proper minimum number of consistent paired-ends . Option \n" "ERROR" >&2;
	usagelongdoc; 
	exit -1; 
    fi
fi
	
# Minimum number of total supporting evidences for novel splice-sites (spanning reads + discordant PE)
if [[ "$minNbTotalNovelSS" == "" ]];
then 
    minNbTotalNovelSS=6;
else
	if [[ ! "$minNbTotalNovelSS" =~ ^[0-9]+$ ]]; 
    then
	log "Please specify a proper minimum supporting evidence for novel splice-sites. Option --total-support-novel-ss\n" "ERROR" >&2;
	usagelongdoc; 
	exit -1; 
    fi
fi

# Minimum number of spanning reads for novel splice-sites
if [[ "$minNbSpanningNovelSS" == "" ]];
then 
    minNbSpanningNovelSS=3;
else
	if [[ ! "$minNbSpanningNovelSS" =~ ^[0-9]+$ ]]; 
    then
	log "Please specify a proper minimum number of spanning reads for novel splice-sites. Option --spanning-reads-novel-ss\n" "ERROR" >&2;
	usagelongdoc; 
	exit -1; 
    fi
fi
	
# Minimum number of consistent paired-ends for novel splice-sites
if [[ "$minNbConsistentPENovelSS" == "" ]];
then 
    minNbConsistentPENovelSS=3;
else
	if [[ ! "$minNbConsistentPENovelSS" =~ ^[0-9]+$ ]]; 
    then
	log "Please specify a proper minimum number of consistent paired-ends for novel splice-sites. Option --consistent-pairs-novel-ss\n" "ERROR" >&2;
	usagelongdoc; 
	exit -1; 
    fi
fi	
	  
# Minimum percentage of staggered reads
if [[ "$minPercStaggered" == "" ]];
then 
    minPercStaggered=0;
else
	if [[ ! "$minPercStaggered" =~ ^[0-9]+$ ]]; 
    then
	log "Please specify a proper minimum percentage of staggered reads. Option --perc-staggered\n" "ERROR" >&2;
	usagelongdoc; 
	exit -1; 
    fi
fi
		      
# Maximum percentage of multimapped spanning reads  
if [[ "$maxPercMultimaps" == "" ]];
then 
    maxPercMultimaps=100;
else
	if [[ ! "$maxPercMultimaps" =~ ^[0-9]+$ ]]; 
    then
	log "Please specify a proper maximum percentage of multimapped spanning reads. Option --perc-multimappings\n" "ERROR" >&2;
	usagelongdoc; 
	exit -1; 
    fi
fi
	
# Maximum percentage of inconsistent paired ends
if [[ "$maxPercInconsistentPE" == "" ]];
then 
    maxPercInconsistentPE=100;
else
	if [[ ! "$maxPercInconsistentPE" =~ ^[0-9]+$ ]]; 
    then
	log "Please specify a proper maximum percentage of inconsistent paired ends. Option --perc-inconsistent-pairs\n" "ERROR" >&2;
	usagelongdoc; 
	exit -1; 
    fi
fi
	
# Similarity filter configuration
if [[ "$similarityConf" == "" ]];
then 
    similarityConf='30+90';
else
	if [[ ! "$similarityConf" =~ ^([0-9]+\+[0-9]+)$ ]]; 
    then
	log "Please specify a proper similarity filter configuration. Option --similarity\n" "ERROR" >&2;
	usagelongdoc; 
	exit -1; 
    fi
fi

# Biotype black list
if [[ "$biotype" == "" ]];
then 
	biotype='pseudogene,polymorphic_pseudogene,IG_C_pseudogene,IG_J_pseudogene,IG_V_pseudogene,TR_J_pseudogene,TR_V_pseudogene';
fi

# Similarity between gene pairs file
if [[ "$simGnPairs" == "" ]]
then
    simGnPairs="NOT_PROVIDED";
elif [ ! -s "$simGnPairs" ]
then 
    log "Your text file containing similarity information between gene pairs in the annotation does not exist. Option --similarity-gene-pairs\n" "ERROR" >&2; 
    usagelongdoc; 
    exit -1; 
fi


# 4. Directories
################
## binaries and scripts
binDir=$rootDir/bin
awkDir=$rootDir/src/awk
bashDir=$rootDir/src/bash

## Output files directories

# 1. Processed annotation
annotDir=$outDir/Annotation

# 2. Mapping phase
mappingPhaseDir=$outDir/MappingPhase
firstMappingDir=$mappingPhaseDir/FirstMapping
secondMappingDir=$mappingPhaseDir/SecondMapping

# 3. Chimera detection phase
chimeraDetPhaseDir=$outDir/ChimeraDetectionPhase
chimSpliceDir=$chimeraDetPhaseDir/ChimSplice
chimPEDir=$chimeraDetPhaseDir/ChimPE

# 4. Gene similarity 
simDir=$outDir/GnSimilarity


# The temporary directory will be exported as an environmental variable since it will 
# be used by every ChimPipe's scripts 
export TMPDIR=$TMPDIR

# 5. Programs/Scripts
#####################
# Bin 
gemtools=$binDir/gemtools-1.7.1-i3/gemtools
gemrnatools=$binDir/gemtools-1.7.1-i3/gem-rna-tools
gtFilterRemove=$binDir/gemtools-1.7.1-i3/gt.filter.remove
gem2sam=$binDir/gemtools-1.7.1-i3/gem-2-sam
gemInfo=$binDir/gemtools-1.7.1-i3/gem-info
pigz=$binDir/pigz


# Bash 
qual=$bashDir/detect.fq.qual.sh
nbMismLen=$bashDir/nbMismatchesReadLen.sh
addXS=$bashDir/sam2cufflinks.sh
infer_library=$bashDir/infer_library_type.sh
ChimSplice=$bashDir/ChimSplice.sh
chimPE=$bashDir/ChimPE.sh
sim=$bashDir/similarity_bt_gnpairs.sh


# Awk 
processAnnot=$awkDir/preprocess_annotation.awk
addMateInfoSam=$awkDir/add_mateInfo_SAM.awk
gff2Gff=$awkDir/gff2gff.awk
ChimIntegrate=$awkDir/integrate_ChimSplice_ChimPE_output.awk
AddSimGnPairs=$awkDir/add_sim_bt_gnPairs.awk
ChimClassifier=$awkDir/ChimClassifier.awk
ChimFilter=$awkDir/ChimFilter.awk


## DISPLAY PIPELINE CONFIGURATION  
##################################
printf "\n"
header="CHIMPIPE CONFIGURATION FOR $lid"
echo $header
eval "for i in {1..${#header}};do printf \"-\";done"
printf "\n\n"
printf "  %-34s %s\n\n" "ChimPipe Version $version"
printf "  %-34s %s\n" "***** MANDATORY ARGUMENTS *****"

if [[ "$bamAsInput" == "FALSE" ]];
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
	printf "  %-34s %s\n" "bam:" "$bamFirstMap"
	printf "  %-34s %s\n" "genome-index:" "$genomeIndex"
	printf "  %-34s %s\n" "annotation:" "$annot"
	printf "  %-34s %s\n\n" "sample-id:" "$lid"

	printf "  %-34s %s\n" "** Reads information **"
	printf "  %-34s %s\n" "seq-library:" "$readDirectionality"
	printf "  %-34s %s\n\n" "max-read-length:" "$maxReadLength"
	
	printf "  %-34s %s\n" "***** MAPPING PHASE *****"
fi

printf "  %-34s %s\n" "** 2nd mapping **"
printf "  %-34s %s\n" "consensus-ss-fm:" "$spliceSitesSM"
printf "  %-34s %s\n" "min-split-size-fm:" "$splitSizeSM"
printf "  %-34s %s\n\n" "refinement-step-size-fm (0:disabled):" "$refinementSM"	
printf "  %-34s %s\n" "***** CHIMERA DETECTION PHASE *****"
printf "  %-34s %s\n" "** Classification **"
printf "  %-34s %s\n\n" "readthrough-max-dist:" "$readthroughMaxDist"
printf "  %-34s %s\n" "** Filters **"
printf "  %-34s %s\n" "total-support:" "$minNbTotal"
printf "  %-34s %s\n" "spanning-reads:" "$minNbSpanning"
printf "  %-34s %s\n" "consistent-pairs:" "$minNbConsistentPE"
printf "  %-34s %s\n" "total-support-novel-ss:" "$minNbTotalNovelSS"
printf "  %-34s %s\n" "spanning-reads-novel-ss:" "$minNbSpanningNovelSS"
printf "  %-34s %s\n" "consistent-pairs-novel-ss:" "$minNbConsistentPENovelSS"
printf "  %-34s %s\n" "perc-staggered (disabled:0):" "$minPercStaggered"
printf "  %-34s %s\n" "perc-multimappings (disabled:100):" "$maxPercMultimaps"
printf "  %-34s %s\n" "perc-inconsistent-pairs (disabled:100):" "$maxPercInconsistentPE"
printf "  %-34s %s\n" "similarity:" "$similarityConf"
printf "  %-34s %s\n\n" "biotype:" "$biotype"
printf "  %-34s %s\n" "** Files **"
printf "  %-34s %s\n\n" "similarity-gene-pairs:" "$simGnPairs"

printf "  %-34s %s\n" "***** GENERAL *****"
printf "  %-34s %s\n" "output-dir:" "$outDir"
printf "  %-34s %s\n" "tmp-dir:" "$TMPDIR"
printf "  %-34s %s\n" "threads:" "$threads"
printf "  %-34s %s\n" "log:" "$logLevel"
printf "  %-34s %s\n\n" "cleanup:" "$cleanup"

###################
## START CHIMPIPE #
###################
header="Executing ChimPipe $version for $lid"
echo $header
eval "for i in {1..${#header}};do printf \"-\";done"
printf "\n\n"
pipelineStart=$(date +%s)


#######################    	
# 0) PRELIMINARY STEPS #
#######################

## 0.1) Make directories
#######################
# Processed annotation
if [[ ! -d $annotDir ]]; then mkdir $annotDir; fi

# Mapping phase
if [[ ! -d $mappingPhaseDir ]]; then mkdir $mappingPhaseDir; fi
if [[ ! -d $secondMappingDir ]]; then mkdir $secondMappingDir; fi

# Chimera detection phase
if [[ ! -d $chimeraDetPhaseDir ]]; then mkdir $chimeraDetPhaseDir; fi
if [[ ! -d $chimSpliceDir ]]; then mkdir $chimSpliceDir; fi
if [[ ! -d $chimPEDir ]]; then mkdir $chimPEDir; fi

## 0.2) Check quality offset if FASTQ as input
#############################################
if [[ "$bamAsInput" == "FALSE" ]];
then
	step="PRELIM"
	log "Determining the offset quality of the reads for ${lid}..." $step
	run "quality=\`$qual $fastq1 | awk '{print \$2}'\`" "$ECHO" 
	log " The read quality is $quality\n"
	log "done\n"
else
	quality="33"
fi

## 0.3) Generate genome file and 
#################################
# retrieve length longest chromosome    	
#####################################

$gemInfo $genomeIndex | awk '$1 ~ /^#/ && ($2 !~ /M/) && ($2 !~ /Mt/) && ($2 !~ /MT/) && ($3 ~ /F/) {split($5,lengths,"@"); chrLengths[$2]=lengths[2];}END{for (chr in chrLengths){print chr, chrLengths[chr]+1;}}' | tr -d \'\" | sort -V -k1,1 > $annotDir/chromosomes_length.txt

longestChrLength=`sort -k2 -n -r $annotDir/chromosomes_length.txt | awk '{print $2}' | head -1`

## 0.4) Process annotation
#########################

awk -f $processAnnot $annot | awk '($1 !~ /M/) && ($1 !~ /Mt/) && ($1 !~ /MT/)' | awk -f $gff2Gff | sort -V -k1,1 -k4,4n -k5,5n | gzip > $annotDir/annotatedExons.gff.gz
    	

####################    	
# 1) MAPPING PHASE #
####################
 	
# 1.1) First mapping step (skipped if BAM as input). Map all the reads to the genome, to the transcriptome and de-novo, using the 
#################################################################################################
# gemtools RNA-Seq pipeline but with max intron size larger than the biggest chromosome, and with 
##################################################################################################  
# a number of mismatches of round(read_length/6) an edit distance of round(read_length/20) 
##########################################################################################
# outputs are: 
##############
# - $outDir/${lid}.map.gz 
# - $outDir/${lid}_raw_chrSorted.bam


if [[ "$bamAsInput" == "FALSE" ]];
then
	gemFirstMap=$firstMappingDir/${lid}_firstMap.map.gz
	bamFirstMap=$firstMappingDir/${lid}_firstMap.bam
	
	if [ ! -s $bamFirstMap ]; 
	then
		if [[ ! -d $firstMappingDir ]] ; then mkdir $firstMappingDir; fi
		firstMapping_FASTQinput		## Call function to run the gemtools rna-pipeline
	else
		printHeader "First mapping BAM file already exists... skipping first mapping step";
	fi
else
	printHeader "BAM file provided as input... skipping first mapping step";
fi


# 1.2) Extract first mapping unmapped reads for a second split-mapping
##########################################################################
# attemp allowing split-mappings in different chromosomes, strands 
######################################################################
# and non genomic order. Produce a FASTQ file with them. 
########################################################
# output is: 
############
# - $outDir/${lid}_reads2remap.fastq

## Comment: samtools view -f 4  ** extract unmapped reads  

reads2remap=$secondMappingDir/${lid}_reads2remap.fastq


if [ ! -s $reads2remap ]; 
then
	step="READS2REMAP"
	startTime=$(date +%s)
	printHeader "Executing extract reads to remap step" 
	run "samtools view -h -f 4 -@ $threads $bamFirstMap | awk -v OFS="'\\\t'" -f $addMateInfoSam | samtools view -@ $threads -bS - | bedtools bamtofastq -i - -fq $reads2remap >> $secondMappingDir/${lid}_reads2remap.log 2>&1" "$ECHO"	

    if [ -s $reads2remap ]; 
    then
        endTime=$(date +%s)
		printHeader "Extracting reads completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
    else	    
        log "Error extracting the reads\n" "ERROR" 
        exit -1
	fi
else
    printHeader "FASTQ file with reads to remap already exists... skipping extracting reads to remap step"
fi

# 1.3) Second split-mapping attemp. Remap the extracted reads allowing reads 
###########################################################################
# to split in different chromosomes, strands and non genomic order.
###################################################################
# output is: 
############
# - $outDir/SecondMapping/${lid}.remapped.map

gemSecondMap=$secondMappingDir/${lid}_secondMap.map

if [ ! -s $gemSecondMap ];
then
	step="SECOND-MAP"
	startTime=$(date +%s)
	printHeader "Executing second split-mapping step"
	log "Remapping reads allowing them to split-map in different chromosomes, strand and non genomic order..." $step	
	run "$gemrnatools split-mapper -I $genomeIndex -i $reads2remap -q 'offset-33' -o ${gemSecondMap%.map} -t 10 -T $threads --min-split-size $splitSizeSM --refinement-step-size $refinementSM --splice-consensus $spliceSitesSM  >> $secondMappingDir/${lid}_secondMap.log 2>&1" "$ECHO"
	
	if [ ! -s $gemSecondMap ]; 
	then
        log "Error in the second mapping\n" "ERROR" 
        exit -1
    fi
    
    endTime=$(date +%s)
	printHeader "Unmapped reads mapping step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
	printHeader "Second mapping GEM file already exists... skipping extracting second mapping step"
fi
	

# 1.4) Infer the sequencing library protocol used (UNSTRANDED, MATE2_SENSE OR MATE1_SENSE) 
########################################################################################
# from a subset with the 1% of the mapped reads. 
#################################################
# Outputs are: 
##############
# - variables $readDirectionality and $stranded

if [[ "$readDirectionality" == "UNKNOWN" ]]; 
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


##############################    	
# 2) CHIMERA DETECTION PHASE #
##############################

# 2.1) run ChimSplice on first and second mappings 
###################################################
# output files: 
################
# - $chimSpliceDir/chimeric_spliceJunctions.txt
# - $chimSpliceDir/normal_spliceJunctions.txt
# - $chimSpliceDir/unannotated_spliceJunctions.txt

paths2ChimSplice=$chimSpliceDir/first_second_mapping_aligment_files_paths_$lid.txt
chimSpliceOut=$chimSpliceDir/chimeric_spliceJunctions.txt

### Execute ChimSplice

if [ ! -s $chimSpliceOut ]; 
then
	step="CHIMSPLICE"
	startTime=$(date +%s)
	printHeader "Executing ChimSplice"
	log "Finding chimeric junctions from split-mappings..." $step
	run "echo $bamFirstMap > $paths2ChimSplice" "$ECHO"
	run "echo $gemSecondMap >> $paths2ChimSplice" "$ECHO"
	run "$ChimSplice $paths2ChimSplice $genomeIndex $annotDir/annotatedExons.gff.gz $readDirectionality $spliceSitesFM $chimSpliceDir 1> $chimSpliceDir/ChimSplice.out 2> $chimSpliceDir/ChimSplice.err" "$ECHO"
	
	if [ ! -s $chimSpliceOut ]; 
	then
        log "Error running ChimSplice\n" "ERROR" 
        exit -1
    fi
    
    endTime=$(date +%s)
	printHeader "ChimSplice step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
	printHeader "Chimeric Junctions file already exists... skipping step"
fi

# 2.2) run ChimPE on first and second mappings 
###################################################
# output files: 
################
# - $chimPEDir/discordant_readPairs.txt
# - $chimPEDir/concordant_contiguousMapped_readPairs.txt
# - $chimPEDir/unannotated_contiguousMapped_readPairs.txt

chimPEOut=$chimPEDir/discordant_readPairs.txt

if [ ! -s $chimPEOut ];
then
	step="CHIMPE"
	startTime=$(date +%s)
	printHeader "Executing ChimPE"
	log "Finding discordant read pairs connecting exons from two different genes..." $step
	run "$chimPE $bamFirstMap $annotDir/annotatedExons.gff.gz $chimSpliceDir/normal_spliceJunctions.txt $readDirectionality $chimPEDir 1> $chimPEDir/ChimPE.out 2> $chimPEDir/ChimPE.err" "$ECHO"
	
	if [ ! -s $chimPEOut ]; 
	then
        log "Error running ChimPE\n" "ERROR" 
        exit -1
    fi
    
    endTime=$(date +%s)
	printHeader "ChimPE step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
	printHeader "Discordant paired-end file already exists... skipping step"
fi

# 2.3) Integrate ChimSplice and ChimPE outputs
##############################################
# output file: 
##############
# - $chimeraDetPhaseDir/chimeric_spliceJunctions_discordantPE.txt

ChimIntegrateOut=$chimeraDetPhaseDir/chimeric_spliceJunctions_discordantPE.txt

if [ ! -s $ChimIntegrateOut ];
then
	step="CHIMINTEGRATE"
	startTime=$(date +%s)
	printHeader "Executing ChimIntegrate"
	log "Integrate ChimSplice and ChimPE output files..." $step
	run "awk -v ChimSplice=<(awk 'NR>1' $chimSpliceOut) -f $ChimIntegrate $chimPEOut > $ChimIntegrateOut" "$ECHO"
	
	if [ ! -s $ChimIntegrateOut ]; 
	then
        log "Error running ChimIntegrate\n" "ERROR" 
        exit -1
    fi
    
    endTime=$(date +%s)
	printHeader "ChimIntegrate step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
	printHeader "ChimIntegrate output file already exists... skipping step"
fi

# 2.4) Compute the gene similarity matrix in case the user does not provide it
###############################################################################
# Output:
# - $outDir/$b2\_gnPair_similarity_matrix.txt
# - $outDir/chimeric_spliceJunctions_candidates_${lid}.txt

if [ ! -e "$simGnPairs" ]
then
 	if [[ ! -d $simDir ]]; then mkdir $simDir; fi
	cd $simDir
    
    printHeader "Executing ChimSimilarity"
    step="CHIMSIM"
    startTime=$(date +%s)
    log "Computing similarity between annotated genes..." $step
    run "$sim $annot $genomeIndex  1> $simDir/sim.out 2> $simDir/sim.err" "$ECHO"
    
    ## Define variable with annotation name:
	b=`basename $annot`
	b2tmp=${b%.gtf}
	b2=${b2tmp%.gff}
    
    ## Output file name:
    simGnPairs=$simDir/$b2.similarity.txt
    
    if [ ! -s $simGnPairs ]; 
	then
        log "Error running ChimSim\n" "ERROR" 
        exit -1
    fi
 	
    endTime=$(date +%s)
    
    cd $outDir
    printHeader "Computing similarity between annotated gene step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
	printHeader "Similarity matrix already exists... skipping step"
fi

### Add information regarding the sequence similarity between connected genes 
# to the chimeric junctions matrix

chimJuncCandidates=$outDir/chimericJunction_candidates_${lid}.txt

if [  ! -e "$chimJuncCandidates" ]
then		
	step="SIM"
	startTime=$(date +%s)
	log "Adding sequence similarity between connected genes information to the chimeric junction candidates..." $step
	run "awk -v fileRef=$simGnPairs -f $AddSimGnPairs $ChimIntegrateOut > $chimJuncCandidates" "$ECHO"
	
	if [ ! -s $chimJuncCandidates ]
	then
	    log "Error adding similarity information\n" "ERROR" 
	    exit -1
	else
		rm $ChimIntegrateOut 
	fi
	
	endTime=$(date +%s)
	printHeader "Add sequence similarity information step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Similarity information already added to chimeric junctions... skipping step"
fi


# 2.5) Classify chimera candidates 
###################################
# - $outDir/chimericJunction_candidates_classified_${lid}.txt

classifiedCandidates=$outDir/chimericJunction_candidates_classified_${lid}.txt

if [  ! -e "$classifiedCandidates" ]
then		
	step="CLASSIFY"
	startTime=$(date +%s)
	log "Classifying chimera candidates..." $step
	run "awk -v readthroughMaxDist=$readthroughMaxDist -f $ChimClassifier $chimJuncCandidates > $classifiedCandidates" "$ECHO"
	
	if [ ! -s $classifiedCandidates ]
	then
	    log "Error classifying candidates\n" "ERROR" 
	    exit -1
	fi
	endTime=$(date +%s)
	printHeader "Classifying chimeric junction candidates step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Chimeric junction candidates already classified... skipping step"
fi

# 2.6) Filter chimera candidates to produce a final set of chimeric junctions
################################################################################
# - $outDir/chimeric_junctions.txt

chimJuncCandidatesFiltered=$outDir/chimericJunction_candidates_filtered_${lid}.txt
chimJunc=$outDir/chimericJunctions_${lid}.txt
chimJuncFiltered=$outDir/chimericJunctions_filtered_${lid}.txt

if [ ! -s "$chimJunc" ]; 
then
	step="CHIMFILTER"
	startTime=$(date +%s)
	printHeader "Executing ChimFilter"
	log "Filtering chimera candidates..." $step	
	run "awk -v minNbTotal=$minNbTotal -v minNbSpanning=$minNbSpanning -v minNbConsistentPE=$minNbConsistentPE -v minNbTotalNovelSS=$minNbTotalNovelSS -v minNbSpanningNovelSS=$minNbSpanningNovelSS -v minNbConsistentPENovelSS=$minNbConsistentPENovelSS -v minPercStaggered=$minPercStaggered -v maxPercMultimaps=$maxPercMultimaps -v maxPercInconsistentPE=$maxPercInconsistentPE -v similarityConf=$similarityConf -v biotype=$biotype -f $ChimFilter $classifiedCandidates > $chimJuncCandidatesFiltered" "$ECHO"
	
	## Make one file with filtered chimeras and another one with the ones that pass the filtering
	awk '(NR==1) || ($3=="0")' $chimJuncCandidatesFiltered > $chimJunc
	awk '(NR==1) || ($3=="1")' $chimJuncCandidatesFiltered > $chimJuncFiltered
	
	if [ ! -s $chimJuncCandidatesFiltered ]; 
	then
		log "Error filtering chimeric junction candidates\n" "ERROR" 
    	exit -1
	else
		rm $classifiedCandidates $chimJuncCandidatesFiltered
	fi 
	
	endTime=$(date +%s)
	printHeader "ChimFilter completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else 
	printHeader "Chimeric junction candidates already filtered... skipping step"
fi

######################
# 3) CLEANUP AND END #
######################

if [[ "$cleanup" == "TRUE" ]]; 
then 
	step="CLEAN-UP"
	startTime=$(date +%s)
	log "Removing intermediate files..." $step
	rm -r $annotDir $firstMappingDir/${lid}_firstMap.map.gz $firstMappingDir/${lid}_firstMap_filtered.map.gz $secondMappingDir/${lid}_reads2remap.fastq $chimeraDetPhaseDir $outDir/chimericJunction_candidates_${lid}.txt 
	log "done\n" 
	endTime=$(date +%s)
	printHeader "Clean up step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else 
	printHeader "No clean up mode... skipping step"
fi	

pipelineEnd=$(date +%s)
printHeader "ChimPipe for $lid completed in $(echo "($pipelineEnd-$pipelineStart)/60" | bc -l | xargs printf "%.2f\n") min "

# disable extglob
shopt -u extglob

exit 0

