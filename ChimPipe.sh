#!/bin/bash

<<authors
*****************************************************************************
	
	ChimPipe.sh
	
	This file is part of the ChimPipe pipeline 

	Copyright (c) 2014 Bernardo Rodríguez-Martín 
					   Emilio Palumbo 
					   Sarah djebali 
	
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

# Function 1. Print stdout basic usage information
###################################################
function usage
{
cat <<help
	
*** ChimPipe version $version ***

Execute ChimPipe (from paired-end RNA-Seq reads to chimeric junctions) on one RNA-Seq dataset (sample).
	
USAGE: $0 -i <fastq_file> -g <genome_index> -a <annotation> [OPTIONS]

** Mandatory arguments:
	
	-i|--input			<INPUT_FILE>	First mate sequencing reads. ChimPipe deals with paired-end data.
                          				Please make sure the second mate file is in the same directory as
                           				the first one, and the files are named according to the same convention.
                           				E.g: the second mate of "reads_1.fastq" should be "reads_2.fastq".
	
	-g|--genome-index		<GEM>		Index for the reference genome in GEM format.
	
	-a|--annotation			<GTF>		Reference genome annotation file in GTF format. The transcriptome
                                       			index has to be in the same directory as the annotation.

** [OPTIONS] can be:

General:
	-e|--sample-id			<STRING>	Sample identifier (the output files will be named according to this id).   
	-o|--output-dir			<PATH>		Output directory. Default current working directory. 
	--tmp-dir			<PATH>		Temporary directory. Default /tmp.
	-t|--threads			<PATH>		Number of threads to use. Default 1.
	--log			<PATH>		Log level [error |warn | info | debug]. Default info.
	--no-cleanup	<FLAG>		Keep intermediate files. 		
	--dry				<FLAG>		Test the pipeline. Writes the command to the standard output.
	-h|--help			<FLAG>		Display partial usage information, only mandatory plus general arguments.
	-f|--full-help			<FLAG>		Display full usage information. 
		
help
}

# Function 2. Print stdout all the other options
#################################################
function usage_long
{
cat <<help
Reads information:
	--max-read-length		<NUMBER>	Maximum read length. This is used to create the de-novo transcriptome and acts as an upper bound. Default 150.
	-l|--seq-library 		<STRING> 	Type of sequencing library [MATE1_SENSE | MATE2_SENSE | UNSTRANDED].
                        				UNSTRANDED for not strand-specific protocol (unstranded data) and the others 
                        				for the different types of strand-specific protocols (stranded data).
Mapping parameters:
	-M|--mism-contiguous-map	<NUMBER>	Maximum number of mismatches for the contiguous mapping steps with the GEM mapper. Default 4?. Not working
	-m|--mism-split-map		<NUMBER>	Maximum number of mismatches for the segmental mapping steps with the GEM rna-mapper. Default 4?.	Not working
	-c|--consensus-splice-sites	<(couple_1)>, ... ,<(couple_s)>	with <couple> := <donor_consensus>+<acceptor_consensus>
                                 			(list of couples of donor/acceptor splice site consensus sequences, default='GT+AG,GC+AG,ATATC+A.,GTATC+AT'
	--min-split-size		<NUMBER>	Minimum split size for the segmental mapping steps. Default 15.
	--stats				<FLAG>		Enable mapping statistics. Default disabled.
	
Chimeric junctions filter:
	--filter-chimeras		<STRING>	Configuration for the filtering module. Quoted string with 4 numbers separated by commas and ended in semicolom, 
							i.e. "1,2,75:50;", where:
											
								1st: minimum number of staggered reads spanning the chimeric junction.
								2nd: minimum number of paired-end reads encompassing the chimeric junction.		
								3rd: maximum similarity between the connected genes.
								4rd: maximum length of the high similar region between the connected genes.
	
							All these conditions have to be fulfilled for a chimeric junction to pass the filter. It is also possible to make 
							complex conditions by setting two different conditions where at least one of them has to be fulfilled. 
							I.e "10,0,0:0;1,1,0:0;". Default "5,0,80:30;1,1,80:30;".	
	--similarity-gene-pairs	<TEXT>			Text file containing similarity information between the gene pairs in the annotation. Needed for the filtering module 
							to discard junctions connecting highly similar genes. If not provided it will be computed inside ChimPipe.
													
help
}

# Function 3. Print stdout a section header for the string variable
#####################################################################
function printHeader {
    string=$1
    echo "`date` ***** $string *****"
}

# Function 4. Print stdout a section header for the string variable
#####################################################################
function printSubHeader {
    string=$1
    echo "`date` * $string *"
}

# Function 5. Print stdout a header for the string variable
############################################################
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

# Function 6. Print stdout a header for the string variable
############################################################
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
                if [[ ! -e $TMPDIR/$annName ]];then
                    log "Copying annotation file to $TMPDIR..." $step
                    run "cp $annot $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "index")
                if [[ ! -e $TMPDIR/`basename $index` ]];then
                    log "Copying genome index file to $TMPDIR..." $step
                    run "cp $index $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "t-index")
                if [[ ! -e $TMPDIR/$annName.gem ]];then
                    log "Copying annotated transcriptome index file to $TMPDIR..." $step
                    run "cp $annot.junctions.gem $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "keys")
                if [[ ! -e $TMPDIR/$annName.junctions.keys ]];then
                    log "Copying annotated transcriptome keys file to $TMPDIR..." $step
                    run "cp $annot.junctions.keys $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            esac
    done
}

# Function 8. Parse user's input
################################
function getoptions {
ARGS=`$getopt -o "i:g:a:e:l:o:t:hfsM:m:c:" -l "input:,genome-index:,annotation:,sample-id:,output-dir:,tmp-dir:,threads:,log:,dry,help,full-help,no-cleanup,mis-contiguous-map:,mism-split-map:,consensus-splice-sites,max-read-length:,seq-library:,max-read-length:,min-split-size:,filter-chimeras:,similarity-gene-pairs:,stats" \
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
        
    -e|--sample-id)
       if [ -n "$2" ];
       then
         lid=$2
       fi
       shift 2;;
      
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

	 -M|--mism-contiguous-map)
       if [ -n "$2" ];
       then
         mism=$2
       fi
       shift 2;;
	
	-m|--mism-split-map)
	    if [ -n "$2" ];
	    then
		mismSplit=$2
	    fi
	    shift 2;;
       
	--min-split-size)
	    if [ -n "$2" ];
	    then
			splitSize=$2
	    fi
	    shift 2;;   
	
	-c|--consensus-splice-sites)
	    if [ -n "$2" ];
	    then
			spliceSites=$2
	    fi
	    shift 2;;
	
	--stats)
	    mapStats="1"  
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
 
    --log)
      	if [ -n $2 ];
      	then
        	logLevel=$2
      	fi
      	shift 2;;
      	 
	-t|--threads)
	    if [ -n $2 ];
	    then
			threads=$2
	    fi
	    shift 2;;
	
	--dry)
	    ECHO="echo "
	    shift;;
	
	-h|--help)
	    usage
	    exit 1
	    shift;;
    
	-f|--full-help)
	    usage
	    usage_long
	    exit 1
	    shift;;
      
	--no-cleanup)
	    cleanup=0;
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
version=V0.8.3

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

# Mandatory arguments
#####################
if [[ ! -e $input ]]; then log "Your input file file does not exist\n" "ERROR" >&2; usage; exit -1; fi
if [[ `basename ${input##*_}` != "1.fastq.gz" ]]; then log "Please check that the name of your FASTQ file ends with \"_1.fastq.gz\"\n" "ERROR" >&2; usage; exit -1; fi
if [[ ! -e $index ]]; then log "Your genome index file does not exist\n" "ERROR" >&2; usage; exit -1; fi
if [[ ! -e $annot ]]; then log "Your annotation file does not exist\n" "ERROR" >&2; usage; exit -1; fi
annName=`basename $annot`


# Optional arguments
#####################

# Maximum read length

if [[ "$maxReadLength" == "" ]]; 
then 
	maxReadLength='150'; 
else
	if [[ ! "$maxReadLength" =~ ^[0-9]+$ ]]; 
	then
		log "Please specify a proper maximum read length value for mapping. Option --max-read-length\n" "ERROR" >&2;
		usage;
		usage_long; 
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
		usage; 
		exit -1;	
	fi
else
	readDirectionality="NOT DEFINED"
fi 	

# Number of mismatches contiguous mapping
if [[ "$mism" == "" ]]; 
then 
	mism='4'; 
else 
	if [[ ! "$mism" =~ ^[0-9]+$ ]]; 
	then
		log "Please specify a proper number of mismatches for contiguous mapping steps. Option -M|--mism-contiguous-map\n" "ERROR" >&2;
		usage; 
		usage_long;
		exit -1; 
	fi
fi

# Consensus splice sites
if [[ "$spliceSites" == "" ]]; 
then 
	spliceSites="GT+AG,GC+AG,ATATC+A.,GTATC+AT"; 
else			
	if [[ ! "$spliceSites" =~ ^([ACGT.]+\+[ACGT.]+,)*([ACGT.]+\+[ACGT.]+)$ ]];
	then
		log "Please specify a proper consensus splice site sequence for the first mapping. Option -c|--consensus-splice-sites\n" "ERROR" >&2;
		usage; 
		usage_long;
		exit -1; 
	fi
fi

# Mapping statistics
if [[ "$mapStats" != "1" ]]; 
then 
	mapStats=0;
	stats="--no-stats"; 
	count="--no-count"; 
fi

# Minimum split size for the segmental mappings
if [[ "$splitSize" == "" ]];
then 
	splitSize='15'; 
else
	if [[ ! "$splitSize" =~ ^[0-9]+$ ]]; 
	then
		log "Please specify a proper minimum split size for the segmental mapping steps. Option --min-split-size\n" "ERROR" >&2;
		usage; 
		usage_long;
		exit -1; 
	fi
fi

# Filtering module configuration	
if [[ "$filterConf" == "" ]]; 
then 			
	filterConf="5,0,80,30;1,1,80,30;";		# Default
else
	if [[ ! "$filterConf" =~ ^([0-9]+,[0-9]+,[0-9]{,3},[0-9]+;){1,2}$ ]]; 
	then
		log "Please check your filtering module configuration. Option --filter-chimeras\n" "ERROR" >&2; 
		usage; 
		usage_long;
		exit -1;
	fi
fi

# Similarity between gene pairs file

if [[ $simGnPairs == "" ]]
then
	simGnPairs="NOT PROVIDED";
elif [[ ! -s $simGnPairs ]]; 
then 
	log "Your text file containing similarity information between gene pairs in the annotation does not exist. Option --similarity-gene-pairs\n" "ERROR" >&2; 
	usage; 
	exit -1; 
fi

# Output directory
if [[ "$outDir" == "" ]]; 
then 
	outDir=${SGE_O_WORKDIR-$PWD};
else
	if [[ ! -e "$outDir" ]]; 
	then
		log "Your output directory does not exist. Option -o|--output-dir\n" "ERROR" >&2;
		usage; 
		exit -1; 
	fi	
fi

# Library id 
if [ ! -n "$lid" ] 
then 			
    lid=myexp		# Default
fi

# Temporary directory
if [[ "$TMPDIR" == "" ]]; 
then 
	TMPDIR='/tmp'; 
else	
	if [[ ! -e "$TMPDIR" ]]; 
	then
		log "Your temporary directory does not exist. Option --tmp-dir\n" "ERROR" >&2;
		usage; 
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
		usage; 
		exit -1; 
	fi
fi

if [[ "$threads" == "1" ]]; 
then 
	hthreads='1'; 
else 
	hthreads=$((threads/2)); 
fi	

# Log level
if [[ "$logLevel" == "" ]]; 
then 
	logLevel='info'; 
else	
	if [[ "$logLevel" != @(error|warn|info|debug) ]];
	then
		log "Please specify a proper log status [error |warn | info | debug]. Option -l|--log\n" "ERROR" >&2;
		usage; 
		exit -1; 
	fi
fi

# Clean up
if [[ "$cleanup" == "" ]]; 
then 
	cleanup='1'; 
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

# Bin 
gemtools=$binDir/gemtools-1.7.1-i3/gemtools
gemrnatools=$binDir/gemtools-1.7.1-i3/gem-rna-tools
gtfilter=$binDir/gemtools-1.7.1-i3/gt.filter
gem2sam=$binDir/gemtools-1.7.1-i3/gem-2-sam
gt_filter=$binDir/gemtools-1.7.1-i3/gt.filter.remove
pigz=$binDir/pigz

# Awk 
unmapped=$awkDir/extract_unmapped.awk
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
header="Pipeline configuration for $lid"
echo $header
eval "for i in {1..${#header}};do printf \"-\";done"
printf "\n\n"
printf "  %-34s %s\n\n" "ChimPipe Version $version"
printf "  %-34s %s\n" "***** Mandatory *****"
printf "  %-34s %s\n" "Input file:" "$input"
printf "  %-34s %s\n" "Reference genome file:" "$index"
printf "  %-34s %s\n\n" "Reference gene annotation file:" "$annot"

printf "  %-34s %s\n" "***** Reads information *****"
printf "  %-34s %s\n" "Sequencing library type:" "$readDirectionality"
printf "  %-34s %s\n" "Maximum read length:" "$maxReadLength"
printf "  %-34s %s\n\n" "Sample identifier:" "$lid"

printf "  %-34s %s\n" "***** Mapping *****"
printf "  %-34s %s\n" "Max number of allowed mismatches contiguous mapping:" "$mism"
printf "  %-34s %s\n" "Max number of allowed mismatches split mapping:" "$mismSplit"
printf "  %-34s %s\n" "Consensus-splice-sites for the first mapping:" "$spliceSites"
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
printf "  %-34s %s\n\n" "Clean up (1:enabled, 0:disabled):" "$cleanup"


## START CHIMPIPE
#################
header="Executing ChimPipe $version for $lid"
echo $header
eval "for i in {1..${#header}};do printf \"-\";done"
printf "\n\n"
pipelineStart=$(date +%s)

# 0) Preliminary steps
######################
step="PRELIM"
log "Determining the offset quality of the reads for ${lid}..." $step
run "quality=\`$qual $input | awk '{print \$2}'\`" "$ECHO" 
log "done\n"
    	
# 1) map all the reads to the genome, to the transcriptome and de-novo, using the gemtools RNA-Seq 
#################################################################################################
#   pipeline but with max intron size larger than the biggest chromosome, and with an edit
##################################################################################################  
#   distance of round(read_length/20) -> get a gem map file and a bam file of "normal" mappings;
################################################################################################
# outputs are: 
##############
# - $outDir/$lid.map.gz 
# - $outDir/$lid\_filtered_cuff.bam

bamFirstMapping=$outDir/${lid}_filtered_cuff.bam

if [ ! -s $bamFirstMapping ]; then
	step="FIRST-MAP"
	startTimeFirstMap=$(date +%s)
	printHeader "Executing first mapping step"    

	# 1.1) Mapping
	################
	if [ ! -s $lid.map.gz ];then
	    step="FIRST-MAP"
	    startTime=$(date +%s)
    
	## Copy needed files to TMPDIR
    	copyToTmp "index,annotation,t-index,keys"

    	log "Running gemtools rna pipeline on ${lid}..." $step
    	run "$gemtools --loglevel $logLevel rna-pipeline -f $input -i $TMPDIR/`basename $index` -a $TMPDIR/$annName -q $quality -n $lid --max-read-length $maxReadLength --max-intron-length 300000000 --junction-consensus $spliceSites -t $threads --no-bam --no-filtered $stats $count" "$ECHO" 
		log "done\n"
    	if [ -s $lid.map.gz ]; then
        	log "Computing md5sum for map file..." $step
        	run "md5sum $lid.map.gz > $lid.map.gz.md5" "$ECHO"
        	log "done\n"
    	else
        	log "Error producing map file\n" "ERROR"
        	exit -1
    	fi
    	endTime=$(date +%s)
    	printSubHeader "Mapping step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
	else
    	printSubHeader "Map file already present...skipping mapping step"
	fi

	# 1.2) Filtering the map file
	##############################
	filteredGem=${lid}_unique_${mism}mism.map.gz

	if [ ! -s $filteredGem ];then
    	step="FIRST-MAP.FILTER"
    	startTime=$(date +%s)
    	printSubHeader "Executing filtering step"
	
    	log "Filtering map file..." $step
    	run "$gt_filter -i $lid.map.gz --max-matches 2 --max-levenshtein-error $mism -t $threads | $pigz -p $threads -c > $filteredGem" "$ECHO"
    	log "done\n" 
    	if [ -s $filteredGem ]; then
    	    log "Computing md5sum for filtered file..." $step
    	    run "md5sum $filteredGem > $filteredGem.md5" "$ECHO"
    	    log "done\n"
    	else
    	    log "Error producing filtered map file\n" "ERROR" 
    	    exit -1
    	fi
    	endTime=$(date +%s)
    	printSubHeader "Filtering step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
	else	
    	printSubHeader "Filtered map file is present...skipping filtering step"
	fi

	# 1.3) Stats from the filtered file
	####################################
	filteredGemStats=${filteredGem%.map.gz}.stats

	if [ $filteredGemStats -ot $filteredGem ] && [ "$mapStats" == "1" ] ;then
    	step="FIRST-MAP.STATS"
    	startTime=$(date +%s)
    	printSubHeader "Executing GEM stats step"
    	log "Producing stats for $filteredGem..." $step
    	run "$gemtools stats -i $filteredGem -t $threads -a -p 2> $filteredGemStats" "$ECHO"
    	log "done\n"
    	if [ -s $filteredGemStats ]; then
        	log "Computing md5sum for stats file..." $step
        	run "md5sum $filteredGemStats > $filteredGemStats.md5" "$ECHO"
        	log "done\n"
    	else
        	log "Error producing GEM stats\n" "ERROR" 
        	exit -1
    	fi
    	endTime=$(date +%s)
    	printSubHeader "GEM stats step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
	else
    	if [ "$mapStats" == "1" ]; then 
    		printSubHeader "GEM stats file is present...skipping GEM stats step"; 
    	fi 
	fi

	# 1.4) Convert to bam and adding the XS field
	###############################################
	filteredBam=${lid}_filtered_cuff.bam

	if [ ! -s $filteredBam ];then
	    step="FIRST-MAP.CONVERT"
	    startTime=$(date +%s)
	    printSubHeader "Executing conversion step"
	
	    ## Copy needed files to TMPDIR
	    copyToTmp "index"	
	    log "Converting $lid to bam..." $step
    	run "$pigz -p $hthreads -dc $filteredGem | $gem2sam -T $hthreads -I $TMPDIR/`basename $index` --expect-paired-end-reads -q offset-$quality -l | sed 's/chrMT/chrM/g' | $addXS $readDirectionality | samtools view -@ $threads -Sb - | samtools sort -@ $threads -m 4G - $TMPDIR/${filteredBam%.bam}" "$ECHO"
    	log "done\n"
    	if [ -s $TMPDIR/$filteredBam ]; then
        	log "Computing md5sum for bam file..." $step
        	run "md5sum $TMPDIR/$filteredBam > $TMPDIR/$filteredBam.md5" "$ECHO"
        	run "cp $TMPDIR/$filteredBam.md5 ." "$ECHO"
        	log "done\n"
        	log "Copying filtered bam file to mapping dir..." $step
        	run "cp $TMPDIR/$filteredBam ." "$ECHO"
        	log "done\n"
    	else
        	log "Error producing the filtered bam file\n" "ERROR" 
        	exit -1
    	fi
    	endTime=$(date +%s)
    	printSubHeader "Conversion step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
	else
    	printSubHeader "Bam file is present...skipping conversion step"
	fi
	endTimeFirstMap=$(date +%s)
	printHeader "First mapping for $lid completed in $(echo "($endTimeFirstMap-$startTimeFirstMap)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "First mapping BAM file is present...skipping first mapping step"
fi

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
	read fraction1 fraction2 other <<<$(bash $infer_library $bamFirstMapping $annot)
	log "done\n"
	log "Fraction of reads explained by 1++,1--,2+-,2-+: $fraction1\n" $step
 	log "Fraction of reads explained by 1+-,1-+,2++,2--: $fraction2\n" $step
	log "Fraction of reads explained by other combinations: $other\n" $step 
	
	# Turn the percentages into integers
	fraction1_int=${fraction1/\.*};
	fraction2_int=${fraction2/\.*};
	other_int=${other/\.*};

	# Infer the sequencing library from the mapping distribution. 
	if [ "$fraction1_int" -ge 80 ]; # MATE1_SENSE protocol
	then 
		readDirectionality="MATE1_SENSE";
		stranded=1;
		echo $readDirectionality;	
	elif [ "$fraction2_int" -ge 80 ];
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
			usage
    		usage_long
			exit -1	
		fi
	else
		log "ChimPipe is not able to determine the library type. Ask your data provider and use the option -l|--seq-library\n" "ERROR" >&2;
		usage
    	usage_long
		exit -1	
	fi
	log "Sequencing library type: $readDirectionality\n" $step 
	log "Strand aware protocol (1: yes, 0: no): $stranded\n" $step 
	endTimeLibrary=$(date +%s)
	printHeader "Sequencing library inference for $lid completed in $(echo "($endTimeLibrary-$startTimeLibrary)/60" | bc -l | xargs printf "%.2f\n") min"
else
	printHeader "Sequencing library type provided by the user...skipping library inference step"
fi

exit 1

# 3) extract the reads that do not map with a number of mismatches lower than 6
###############################################################################
#    from the gem file: outDirget a fastq file 
########################################
# output is: 
############
# - $outDir/$lid.unmapped.fastq

unmappedReads=$outDir/${lid}.unmapped.fastq

if [ ! -s $unmappedReads ]; then
	step="UNMAP"
	startTime=$(date +%s)
	printHeader "Executing unmapped reads step" 
	log "Extracting the reads that do not map with a number of mismatches lower than 6..." $step
	run "zcat $lid.map.gz | awk -f $unmapped | $gtfilter -t $threads --output-format 'FASTA' > $outDir/$lid.unmapped.fastq" "$ECHO" 	
	log "done\n" 
    if [ -s $unmappedReads ]; then
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

# 4) map the unmapped reads with the rna mapper binary (!!!parameters to think!!!): get a gem map file
######################################################################################################
#   of "atypical" mappings
########################
# output is: 
############
# - $outDir/SecondMapping/$lid.unmapped_rna-mapped.map

gemSecondMapping=$outDir/SecondMapping/${lid}.unmapped_rna-mapped.map

if [ ! -s $gemSecondMapping ]; then
	step="SECOND-MAP"
	startTime=$(date +%s)
	printHeader "Executing second mapping step"
	log "Mapping the unmapped reads with the rna mapper..." $step	
	run "$gemrnatools split-mapper -I $index -i $outDir/$lid.unmapped.fastq -q 'offset-$quality' -o $outDir/SecondMapping/$lid.unmapped_rna-mapped -t 10 -T $threads -c 'GT+AG' --min-split-size $splitSize > $outDir/SecondMapping/$lid.gem-rna-mapper.out" "$ECHO"
	log "done\n" 
	if [ -s $gemSecondMapping ]; then
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

# 5) extract the reads mapping both uniquely and in 2 blocks from the bam file of "normal" mappings and convert in gff.gz
#########################################################################################################################
# output is: 
############
# - $outDir/FromFirstBam/$lid\_filtered_cuff_2blocks.gff.gz

gffFromBam=$outDir/FromFirstBam/${lid}_filtered_cuff_2blocks.gff.gz

if [ ! -s $gffFromBam ]; then
	step="FIRST-CONVERT"
	startTime=$(date +%s)
	printHeader "Executing conversion of the bam into gff step"
	log "Generating a ".gff.gz" file from the normal mappings containing the reads split-mapping both uniquely and in 2 blocks..." $step
	bedtools bamtobed -i $outDir/$lid\_filtered_cuff.bam -bed12 | awk '$10==2' | awk -v rev='1' -f $bed2bedPE | awk -v readDirectionality=$readDirectionality  -f $bedPECorrectStrand | awk -f $bedPE2gff | awk -f $gff2Gff | gzip > $outDir/FromFirstBam/$lid\_filtered_cuff_2blocks.gff.gz
	log "done\n"
	if [ -s $gffFromBam ]; then
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

# 6) extract the reads mapping both uniquely and in 2 blocks from the map file of "atypical" mappings and convert in gff.gz
#########################################################################################################################
# output is: 
############
# - $outDir/FromSecondMapping/${lid}.unmapped_rna-mapped.gff.gz

gffFromMap=$outDir/FromSecondMapping/${lid}.unmapped_rna-mapped.gff.gz

if [ ! -s $gffFromMap ]; then
	step="SECOND-CONVERT"
	startTime=$(date +%s)
	printHeader "Executing conversion of the gem into gff step"
	log "Generating a ".gff.gz" file from the atypical mappings containing the reads split-mapping both uniquely and in 2 blocks..." $step	
	run "awk -v readDirectionality=$readDirectionality -f $gemCorrectStrand $outDir/SecondMapping/$lid.unmapped_rna-mapped.map | awk -v rev="0" -f $gemToGff | awk -f $gff2Gff | gzip > $outDir/FromSecondMapping/${lid}.unmapped_rna-mapped.gff.gz" "$ECHO"
	log "done\n" 
	if [ -s $gffFromMap ]; then
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


# 7) put the path to the "normal" and "atypical" gff.gz files in a same txt file for chimsplice 
#############################################################################################

run "echo $outDir/FromFirstBam/$lid\_filtered_cuff_2blocks.gff.gz > $outDir/split_mapping_file_sample_$lid.txt" "$ECHO"
run "echo $outDir/FromSecondMapping/$lid.unmapped_rna-mapped.gff.gz >> $outDir/split_mapping_file_sample_$lid.txt" "$ECHO"

# 8) run chimsplice on 4) and 5)
###############################
# - $outDir/Chimsplice/chimeric_junctions_report_$lid.txt
# - $outDir/Chimsplice/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt
# - $outDir/Chimsplice/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn_morethan10staggered.txt

exonConnections1=$outDir/Chimsplice/exonA_exonB_with_splitmapping_part1overA_part2overB_readlist_sm1list_sm2list_staggeredlist_totalist_$lid\_filtered_cuff_2blocks.gff.txt.gz 
exonConnections2=$outDir/Chimsplice/exonA_exonB_with_splitmapping_part1overA_part2overB_readlist_sm1list_sm2list_staggeredlist_totalist_$lid.unmapped_rna-mapped.gff.txt.gz

printHeader "Executing Chimsplice step"
if [ ! -s $exonConnections1 ] || [ ! -s $exonConnections2 ]; then
	step="CHIMSPLICE"
	startTime=$(date +%s)
	log "Finding exon to exon connections from the ".gff.gz" files containing the "normal" and "atypical" mappings..." $step
	run "$chim1 $outDir/split_mapping_file_sample_$lid.txt $annot $outDir/Chimsplice $stranded" "$ECHO"
	log "done\n" 
	if [ ! -s $exonConnections1 ] || [ ! -s $exonConnections2 ]; then
        log "Error running chimsplice\n" "ERROR" 
        exit -1
    fi
    endTime=$(date +%s)
	printHeader "Find exon to exon connections step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
	printHeader "Exon to exon connections file present... skipping step"
fi

chimJunctions=$outDir/Chimsplice/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt
if [ ! -s $chimJunctions ]; then
	step="CHIMSPLICE"
	startTime=$(date +%s)
	log "Finding chimeric junctions from exon to exon connections..." $step
	run "$chim2 $outDir/split_mapping_file_sample_$lid.txt $index $annot $outDir/Chimsplice $stranded $spliceSites > $outDir/Chimsplice/chimeric_junctions_report_$lid.txt 2> $outDir/Chimsplice/find_chimeric_junctions_from_exon_to_exon_connections_$lid.err" "$ECHO"
	log "done\n" 
	if [ ! -s $chimJunctions ]; then
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
if [ ! -s $PEinfo ];then
	step="PAIRED-END"
	startTime=$(date +%s)
	log "Finding gene to gene connections supported by paired-end mappings from the ".bam" containing reads mapping in a unique and continuous way..." $step
	run "$findGeneConnections $outDir/$lid\_filtered_cuff.bam $annot $outDir/PE $readDirectionality" "$ECHO"
	log "done\n" 
	if [ ! -s $PEinfo ]; then
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

if [ ! -s $chimJunctionsPE ];then
	step="PAIRED-END"
	startTime=$(date +%s)
	log "Adding PE information to the matrix containing chimeric junction candidates..." $step
	run "awk -v fileRef=$PEinfo -f $addPEinfo $chimJunctions 1> $outDir/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_ss1_ss2_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_PEinfo_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt
" "$ECHO"
	log "done\n" 
	if [ ! -s $chimJunctionsPE ]; then
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

if [  ! -s "$chimJunctionsSim" ]; then		
	if [ "$simGnPairs" != "NOT PROVIDED" ];then
		step="SIM"
		startTime=$(date +%s)
		log "Adding sequence similarity between connected genes information to the chimeric junction matrix..." $step
		run "awk -v fileRef=$simGnPairs -f $AddSimGnPairs $chimJunctionsPE 1> $chimJunctionsSim" "$ECHO"
		log "done\n" 
		if [ ! -s $chimJunctionsSim ]; then
			log "Error adding similarity information\n" "ERROR" 
	    	exit -1
		fi
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

if [ ! -s "$chimJunctionsCandidates" ]; then		
	step="HEADER"
	log "Adding a header to the matrix containing the chimeric junction candidates..." $step
	if [ -s "$chimJunctionsSim" ]; then		
		run "awk 'BEGIN{print \"juncId\", \"nbstag\", \"nbtotal\", \"maxbeg\", \"maxEnd\", \"samechr\", \"samestr\", \"dist\", \"ss1\", \"ss2\", \"gnlist1\", \"gnlist2\", \"gnname1\", \"gnname2\", \"bt1\", \"bt2\", \"PEsupport\", \"maxSim\", \"maxLgal\";}{print \$0;}' $chimJunctionsSim 1> $chimJunctionsCandidates" "$ECHO"
		
		log "done\n"
	else
		if [ -s "$chimJunctionsPE" ]; then
			run "awk 'BEGIN{print \"juncId\", \"nbstag\", \"nbtotal\", \"maxbeg\", \"maxEnd\", \"samechr\", \"samestr\", \"dist\", \"ss1\", \"ss2\", \"gnlist1\", \"gnlist2\", \"gnname1\", \"gnname2\", \"bt1\", \"bt2\", \"PEsupport\";}{print \$0;}' $chimJunctionsPE 1> $chimJunctionsCandidates" "$ECHO"
			
			log "done\n" 	
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

if [ ! -s "$chimJunctions" ]; then
	step="FILTERING MODULE"
	startTime=$(date +%s)
	log "Filtering out chimera candidates to produce a final set of chimeric junctions..." $step
	awk -v filterConf=$filterConf -f $juncFilter $chimJunctionsCandidates 1> $chimJunctions
	log "done\n" 
	if [ ! -s $chimJunction ]; then
		log "Error filtering chimeric junction candidates\n" "ERROR" 
    	exit -1
	fi
	
	endTime=$(date +%s)
	printHeader "Filtering module step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else 
	printHeader "Chimeric junction candidates already filtered... skipping step"
fi

# 12) Clean up
###############

if [[ "$cleanup" == "1" ]]; then 
	step="CLEAN-UP"
	startTime=$(date +%s)
	log "Removing intermediate files..." $step
	rm -r  $outDir/FromFirstBam $outDir/FromSecondMapping $outDir/Chimsplice
	rm $outDir/split_mapping_file_sample_$lid.txt $chimJunctionsSim $chimJunctionsCandidates $chimJunctionsPE
	log "done\n" 
	endTime=$(date +%s)
	printHeader "Clean up step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else 
	printHeader "No clean up mode... skipping step"
fi	


# 13) END
#########
pipelineEnd=$(date +%s)
printHeader "Chimera Mapping pipeline for $lid completed in $(echo "($pipelineEnd-$pipelineStart)/60" | bc -l | xargs printf "%.2f\n") min "

# disable extglob
shopt -u extglob

exit 0

