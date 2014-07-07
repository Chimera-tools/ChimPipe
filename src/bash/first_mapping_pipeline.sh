#!/bin/bash


function usage {
    echo "Usage: $0 -i <fastq_file> -s <sex> [OPTION]..."
    echo "Execute the Blueprint pipeline (only bam generation) on one sample."
    echo ""
    printf "\t-i\tinput file\n"
    printf "\t-g\tindex for the reference genome\n"
    printf "\t-a\treference annotation\n"
    printf "\t-q\tspecify the quality offset of the dataset [33 | 64 | ignore]\n"
    echo ""
    echo "Options:"
    printf "\t-M\tMax number of mismatches. Default \"4\"\n"
    printf "\t-m\tFlag to enable mapping stats. Default disabled.\n"
    printf "\t-L \tMax read length. This is used to create the de-novo transcriptome and acts as an upper bound. Default \"150\"\n"
    printf "\t-s\tFlag to specify whether the sample has strandness information for the reads. Default \"false\"\n" 
    printf "\t-d\tdirectionality of the reads (MATE1_SENSE, MATE2_SENSE, NONE). Default \"NONE\".\n"
    printf "\t-e\tExperiment identifier (the output files will be named according this id). If not specified, the name is inferred from the input files.\n"
    printf "\t-l\tLog level (error, warn, info, debug). Default \"info\".\n"
    printf "\t-T\tTest the pipeline. Writes the command to the standard output.\n"
    printf "\t-t\tNumber of threads to use. Default 1.\n"
    exit 0
}

function printHeader {
    string=$1
    echo "`date` * $string *"
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

function copyToTmp {
    IFS=',' read -ra files <<< "$1"
    for i in ${files[@]};do
        case $i in
            "annotation")
                if [[ ! -e $TMPDIR/$annName ]];then
                    log "Copying annotation file to $TMPDIR..." $step
                    run "cp $annotation $TMPDIR" "$ECHO"
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
                    run "cp $annotation.gem $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "keys")
                if [[ ! -e $TMPDIR/$annName.junctions.keys ]];then
                    log "Copying annotated transcriptome keys file to $TMPDIR..." $step
                    run "cp $annotation.junctions.keys $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            esac
    done
}

## Parsing arguments
#

while getopts ":i:g:a:q:M:mL:sd:e:l:t:Th" opt; do
  case $opt in
    i)
      input="$OPTARG"
      ;;
    g)
      index="`dirname $OPTARG | xargs readlink -e`/`basename $OPTARG`"
      ;;
    a)
      annotation="`dirname $OPTARG | xargs readlink -e`/`basename $OPTARG`"
      ;;
    q)
      quality=$OPTARG
      ;;
    M)
      mism=$OPTARG
      ;;
    m)
      mapStats='1' 
      ;;
    L)
      maxReadLength=$OPTARG 
      ;; 
    s)
      stranded=1
      ;;
    d)
      readStrand=$OPTARG
      ;;
    e)
      sample=$OPTARG
      ;;  
    l)
      loglevel=$OPTARG
      ;;
    t)
      threads=$OPTARG
      ;;
    T)
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


# SETTING UP THE ENVIRONMENT
############################

# = Directories = #
# Environmental variables 
# rootDir - path to the root folder of ChimPipe pipeline. 
# TMPDIR  - temporary directory
# They are environmental variable defined and exported in the main script

binDir=$rootDir/bin
bashDir=$rootDir/src/bash

# = Binaries = #
addXS=$bashDir/sam2cufflinks.sh
gemtools=$binDir/gemtools-1.7.1-i3/bin/gemtools
gem2sam=$binDir/gemtools-1.7.1-i3/bin/gem-2-sam
samtools=$binDir/samtools-0.1.19/samtools
gt_filter=$binDir/miscellaneous/gt.filter.remove
pigz=$binDir/miscellaneous/pigz

# = Variables and input files = #
if [[ "$input" == "" ]]; then log "Please specify the input file\n" "ERROR" >&2; exit -1; fi
if [[ "$index" == "" ]]; then log "Please specify the genome index file\n" "ERROR" >&2; exit -1; fi
if [[ "$annotation" == "" ]]; then log "Please specify the annotation file\n" "ERROR" >&2; exit -1; fi
if [[ "$quality" == "" ]]; then log "Please specify the quality\n" "ERROR" >&2; exit -1; fi
if [[ "$stranded" == "" ]]; then stranded="0"; fi
if [[ "$readStrand" == "" ]]; then readStrand="NONE"; fi
if [[ "$sample" == "" ]]; then basename=$(basename $input); sample=${basename%_1*}; fi
if [[ "$loglevel" == "" ]]; then loglevel="info"; fi
if [[ "$mism" == "" ]]; then mism="4"; fi
if [[ "$mapStats" != "1" ]]; then stats="--no-stats"; count="--no-count"; fi
if [[ "$maxReadLength" == "" ]]; then maxReadLength="150"; fi
annName=`basename $annotation`
if [[ "$threads" == "" ]]; then threads='1'; fi
if [[ "$threads" == "1" ]]; then hthreads='1'; else hthreads=$((threads/2)); fi	
#threads=${NSLOTS-1}


## START
########
pipelineStart=$(date +%s)

## Mapping
#

if [ ! -e $sample.map.gz ];then
    step="FIRST-MAP"
    startTime=$(date +%s)
    
    ## Copy needed files to TMPDIR
    copyToTmp "index,annotation,t-index,keys"

    log "Running gemtools rna pipeline on ${sample}..." $step
    run "$gemtools --loglevel $loglevel rna-pipeline -f $input -i $TMPDIR/`basename $index` -a $TMPDIR/$annName -q $quality -n $sample --max-read-length $maxReadLength --max-intron-length 300000000 -t $threads --no-bam --no-filtered $stats $count" "$ECHO" 
	log "done\n"
    if [ -e $sample.map.gz ]; then
        log "Computing md5sum for map file..." $step
        run "md5sum $sample.map.gz > $sample.map.gz.md5" "$ECHO"
        log "done\n"
    else
        log "Error producing map file\n" "ERROR"
        exit -1
    fi
    endTime=$(date +%s)
    printHeader "Mapping step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Map file already present...skipping mapping step"
fi


## Filtering the map file
##
filteredGem=${sample}_unique_${mism}mism.map.gz

if [ ! -e $filteredGem ];then
    step="FIRST-MAP.FILTER"
    startTime=$(date +%s)
    printHeader "Executing filtering step"

    log "Filtering map file..." $step
    run "$gt_filter -i $sample.map.gz --max-matches 2 --max-levenshtein-error $mism -t $threads | $pigz -p $threads -c > $filteredGem" "$ECHO"
    log "done\n" 
    if [ -e $filteredGem ]; then
        log "Computing md5sum for filtered file..." $step
        run "md5sum $filteredGem > $filteredGem.md5" "$ECHO"
        log "done\n"
    else
        log "Error producing filtered map file\n" "ERROR" 
        exit -1
    fi
    endTime=$(date +%s)
    printHeader "Filtering step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Filtered map file is present...skipping filtering step"
fi

## Stats from the filtered file
##
filteredGemStats=${filteredGem%.map.gz}.stats

if [ $filteredGemStats -ot $filteredGem ] && [ "$mapStats" == "1" ] ;then
    step="FIRST-MAP.STATS"
    startTime=$(date +%s)
    printHeader "Executing GEM stats step"
    log "Producing stats for $filteredGem..." $step
    run "$gemtools stats -i $filteredGem -t $threads -a -p 2> $filteredGemStats" "$ECHO"
    log "done\n"
    if [ -e $filteredGemStats ]; then
        log "Computing md5sum for stats file..." $step
        run "md5sum $filteredGemStats > $filteredGemStats.md5" "$ECHO"
        log "done\n"
    else
        log "Error producing GEM stats\n" "ERROR" 
        exit -1
    fi
    endTime=$(date +%s)
    printHeader "GEM stats step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    if [ "$mapStats" == "1" ]; then printHeader "GEM stats file is present...skipping GEM stats step"; fi 
fi

## Convert to bam and adding the XS field
##
filteredBam=${sample}_filtered_cuff.bam

if [ ! -e $filteredBam ];then
    step="FIRST-MAP.CONVERT"
    startTime=$(date +%s)
    printHeader "Executing conversion step"

    ## Copy needed files to TMPDIR
    copyToTmp "index"

    log "Converting $sample to bam..." $step
    run "$pigz -p $hthreads -dc $filteredGem | $gem2sam -T $hthreads -I $TMPDIR/`basename $index` --expect-paired-end-reads -q offset-$quality -l | sed 's/chrMT/chrM/g' | $addXS $readStrand | $samtools view -@ $threads -Sb - | $samtools sort -@ $threads -m 4G - $TMPDIR/${filteredBam%.bam}" "$ECHO"
    log "done\n"
    if [ -e $TMPDIR/$filteredBam ]; then
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
    printHeader "Conversion step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Bam file is present...skipping conversion step"
fi
exit 0
