#!/bin/bash
#
#Forcing bash shell
#$ -S /bin/bash
#$ -cwd
#$ -o $JOB_NAME.out
#$ -e $JOB_NAME.err
# will exit if there is an error or in a pipe


set -e -o pipefail
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
    printf "\t-m\tMax number of mismatches. Default \"4\"\n"
    printf "\t-M \tMax read length. This is used to create the de-novo transcriptome and acts as an upper bound. Default \"150\"\n"
    printf "\t-s\tflag to specify whether the sample has strandness information for the reads. Default \"false\"\n"
    printf "\t-d\tdirectionality of the reads (MATE1_SENSE, MATE2_SENSE, NONE). Default \"NONE\".\n"
    printf "\t-l\tLog level (error, warn, info, debug). Default \"info\".\n"
    printf "\t-t\tTest the pipeline. Writes the command to the standard output.\n"
    exit 0
}

function printHeader {
    string=$1
    echo "`date` *** $string ***"
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

function copyToTmp {
    IFS=',' read -ra files <<< "$1"
    for i in ${files[@]};do
        case $i in
            "annotation")
                if [ ! -e $TMPDIR/$annName ];then
                    log "Copying annotation file to $TMPDIR..." $step
                    run "cp $annotation $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "index")
                if [ ! -e $TMPDIR/`basename $index` ];then
                    log "Copying genome index file to $TMPDIR..." $step
                    run "cp $index $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "t-index")
                if [ ! -e $TMPDIR/$annName.gem ];then
                    log "Copying annotated transcriptome index file to $TMPDIR..." $step
                    run "cp $annotation.gem $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "keys")
                if [ ! -e $TMPDIR/$annName.junctions.keys ];then
                    log "Copying annotated transcriptome keys file to $TMPDIR..." $step
                    run "cp $annotation.junctions.keys $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "fastq")
                if [ ! -e $TMPDIR/`basename $input` ];then
                    log "Copying fastq files to $TMPDIR..." $step
                    run "cp $input $TMPDIR" "$ECHO"
                    run "cp ${input/_1.fastq/_2.fastq} $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "map.gz")
                if [ ! -e $TMPDIR/$sample.map.gz ]; then
                    log "Copying map file to $TMPDIR..." $step
                    run "cp $sample.map.gz $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "bam")
                if [ ! -e $TMPDIR/$sample.bam ]; then
                    log "Copying bam file to $TMPDIR..." $step
                    run "cp $sample.bam $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "filtered-bam")
                if [ ! -e $TMPDIR/$filteredBam ]; then
                    log "Copying filtered bam file to $TMPDIR..." $step
                    run "cp $filteredBam $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "filtered-bai")
                if [ ! -e $TMPDIR/$filteredBam.bai ];then
                    log "Copying index for filtered bam file to $TMPDIR..." $step
                    run "cp $filteredBam.bai $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "flux-profile")
                if [ ! -e $TMPDIR/$sample.profile ];then
                    log "Copying profile file to $TMPDIR..." $step
                    run "cp $quantDir/$sample/$sample.profile $TMPDIR" "$ECHO"
                    log "done\n"
                fi
                ;;
            "flux-gtf")
               if [ ! -e $TMPDIR/$sample.gtf ];then
                   log "Copying Flux file to $TMPDIR..." $step
                   run "cp $quantDir/$sample/$sample.gtf $TMPDIR" "$ECHO"
                   log "done\n"
               fi
            esac
    done
}

## Parsing arguments
#

while getopts ":i:m:M:g:q:a:std:l:ph" opt; do
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
    m)
      mism=$OPTARG
      ;;
    M)
      maxReadLength=$OPTARG 
      ;; 
    s)
      stranded=1
      ;;
    t)
      ECHO="echo "
      ;;
    d)
      readStrand=$OPTARG
      ;;
    l)
      loglevel=$OPTARG
      ;;
    p)
      profile=1
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


## Setting up the environment
#
BASEDIR=`dirname ${SGE_O_WORKDIR-$PWD}`
BINDIR=~sdjebali/ENCODE_AWG/Analyses/Mouse_Human/Chimeras/bin/
export PATH=/software/rg/el6.3/gemtools/bin:/software/rg/el6.3/flux-capacitor-1.2.4-SNAPSHOT/bin:$HOME/bin:$PATH

## Setting variables and input files
##
if [[ $input == "" ]];then
    log "Please specify the input file\n" "ERROR" >&2
    exit -1
fi

if [[ $index == "" ]];then
    log "Please specify the genome index file\n" "ERROR" >&2
    exit -1
fi

if [[ $annotation == "" ]];then
    log "Please specify the annotation file\n" "ERROR" >&2
    exit -1
fi

if [[ $quality == "" ]];then
    log "Please specify the quality\n" "ERROR" >&2
    exit -1
fi

if [[ $stranded == "" ]];then
   stranded="0"
fi

if [[ $readStrand == "" ]];then
    readStrand="NONE"
fi

if [[ $loglevel == "" ]];then
    loglevel="info"
fi

if [[ $mism == "" ]];then
    mism="4"
fi

if [[ $maxReadLength == "" ]];then
   maxReadLength="150"
fi



## Test

#echo $input
#echo $index
#echo $annotation
#echo $quality
#echo $mism
#echo $maxReadLength
#echo $stranded
#echo $ECHO
#echo $readStrand
#echo $loglevel
#echo $profile


basename=$(basename $input)
sample=${basename%_1*}
threads=${NSLOTS-1}

annName=`basename $annotation`

## Binaries
#
gem2sam="$BINDIR/gem-2-sam"
samtools="$BINDIR/samtools"
addXS="$BINDIR/sam2cufflinks.sh"
trToGn="$BINDIR/TrtoGn_RPKM.sh"
trToEx="$BINDIR/TrtoEx_RPKM.sh"
bamToContigs="$BINDIR/bamToContigs.sh"
gt_quality="$BINDIR/gt.quality"
gt_filter="$BINDIR/gt.filter.remove"
gt_stats="$BINDIR/gt.stats"
pigz="$BINDIR/pigz"
BAMFLAG="/users/rg/dmitri/bamflag/trunk/bamflag"
makecontig="$BINDIR/contigsNew.py"

hthreads=$((threads/2))
if [[ $hthreads == 0 ]];then
    hthreads=1
fi

## START
#
printHeader "Starting Blueprint pipeline (only bam file generation) for $sample"
pipelineStart=$(date +%s)

## Mapping
#
#if [ ! -e $sample.map.gz ];then
if [ ! -e $sample.stats.all.json ];then
    step="MAP"
    startTime=$(date +%s)
    printHeader "Executing mapping step"

    ## Activate the python virtualenv
    #run ". $BASEDIR/venv/bin/activate" "$ECHO"

    ## Copy needed files to TMPDIR
    copyToTmp "index,annotation,t-index,keys"

    log "Running gemtools rna pipeline on ${sample}" $step
    run "gemtools --loglevel $loglevel rna-pipeline -f $input -i $TMPDIR/`basename $index` -a $TMPDIR/$annName -q $quality --max-read-length $maxReadLength --max-intron-length 300000000 -t $threads --no-bam " "$ECHO" 
    #gemtools --loglevel $loglevel rna-pipeline -f $TMPDIR/$basename -i $TMPDIR/`basename $index` -a $TMPDIR/$annName -t $threads -o $TMPDIR --no-bam
    #gemtools --loglevel $loglevel rna-pipeline -f $TMPDIR/$basename -i $TMPDIR/`basename $index` -a $TMPDIR/$annName -m 150 -t $threads -o $TMPDIR --no-sam

    if [ -f $TMPDIR/${sample}.map.gz ]; then
        log "Computing md5sum for map file..." $step
        run "md5sum $TMPDIR/$sample.map.gz > $TMPDIR/$sample.map.gz.md5" "$ECHO"
        run "cp $TMPDIR/$sample.map.gz.md5 ." "$ECHO"
        log "done\n"
        log "Copying map file..." $step
        run "cp $TMPDIR/${sample}.map.gz ." "$ECHO"
        log "done\n"
    #else
    #    log "Error producing map file" "ERROR" >&2
    #    exit -1
    fi
    endTime=$(date +%s)
    printHeader "Mapping step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Map file already present...skipping mapping step"
fi


## Converting to bam
##

# if [ ! -e $sample.bam ];then
#     step="CONVERT"
#     startTime=$(date +%s)
#     printHeader "Executing conversion step"

#     ## Copy needed files to TMPDIR
#     copyToTmp "index"

#     log "Converting ${sample} to bam\n" $step

#     run "pigz -p $hthreads -dc $sample.map.gz | $gem2sam -T $hthreads -I $TMPDIR/`basename $index` --expect-paired-end-reads -q offset-33 | $samtools view -@ $hthreads -Sb - | $samtools sort -@ $hthreads -m `echo $((4<<30))` - $sample" "$ECHO"
#     #pigz -p $threads -dc $TMPDIR/$sample.map.gz | $gem2sam -T $hthreads -I $TMPDIR/`basename $index` --expect-paired-end-reads -q offset-33 | $samtools view -Sb - | $samtools sort -m $((8<<30)) - $TMPDIR/$sample
#     if [ -f $TMPDIR/${sample}.bam ]; then
#         log "Computing md5sum for bam file..." $step
#         run "md5sum $TMPDIR/$sample.bam > $TMPDIR/$sample.bam.md5" "$ECHO"
#         run "cp $TMPDIR/$sample.bam.md5 ." "$ECHO"
#         log "done\n"

#         log "Copying bam file to mapping dir..." $step
#         run "cp $TMPDIR/${sample}.bam ." "$ECHO"
#         log "done\n"
#     #else
#     #    log "Error producing bam file" "ERROR" >&2
#     #    exit -1
#     fi
#     endTime=$(date +%s)
#     printHeader "Conversion step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
# else
#     printHeader "Bam file already present...skipping conversion step"
# fi

## Filtering the map file
##
filteredGem=${sample}_unique_${mism}mism.map.gz

if [ ! -e $filteredGem ];then
    step="FILTER"
    startTime=$(date +%s)
    printHeader "Executing filtering step"

    log "Filtering map file..." $step
    run "$gt_filter -i $sample.map.gz --max-matches 2 -t $threads --max-levenshtein-error $mism -t $threads | $pigz -p $threads -c > $filteredGem" "$ECHO"
    log "done\n" $step
    if [ -f $filteredGem ]; then
        log "Computing md5sum for filtered file..." $step
        run "md5sum $filteredGem > $filteredGem.md5" "$ECHO"
        log "done\n"
    # else
    #    log "Error producing filtered map file" "ERROR" >&2
    #    exit -1
    fi
    endTime=$(date +%s)
    printHeader "Filtering step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Filtered map file is present...skipping fltering step"
fi

## Stats from the filtered file
##
filteredGemStats=${filteredGem%.map.gz}.stats

if [ $filteredGemStats -ot $filteredGem ];then
    step="GEM-STATS"
    startTime=$(date +%s)
    printHeader "Executing GEM stats step"

    ## Copy needed files to TMPDIR
    # copyToTmp "index"

    log "Producing stats for $filteredGem..." $step
    run "$gt_stats -i $filteredGem -t $threads -a -p 2> $filteredGemStats" "$ECHO"
    log "done\n" $step
    if [ -f $filteredGemStats ]; then
        log "Computing md5sum for stats file..." $step
        run "md5sum $filteredGemStats > $filteredGemStats.md5" "$ECHO"
        log "done\n"
    # else
    #    log "Error producing GEM stats" "ERROR" >&2
    #    exit -1
    fi
    endTime=$(date +%s)
    printHeader "GEM stats step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "GEM stats file is present...skipping GEM stats step"
fi

## Convert to bam and adding the XS field
##
filteredBam=${sample}_filtered_cuff.bam

if [ ! -e $filteredBam ];then
    step="CONVERT"
    startTime=$(date +%s)
    printHeader "Executing conversion step"

    ## Copy needed files to TMPDIR
    copyToTmp "index"

    log "Converting  $sample to bam..." $step
    run "$pigz -p $hthreads -dc $filteredGem | $gem2sam -T $hthreads -I $TMPDIR/`basename $index` --expect-paired-end-reads -q offset-$quality -l | sed 's/chrMT/chrM/g' | $addXS $readStrand | $samtools view -@ $threads -Sb - | $samtools sort -@ $threads -m 4G - $TMPDIR/${filteredBam%.bam}" "$ECHO"
    log "done\n" $step
    if [ -f $TMPDIR/$filteredBam ]; then
        log "Computing md5sum for filtered file..." $step
        run "md5sum $TMPDIR/$filteredBam > $TMPDIR/$filteredBam.md5" "$ECHO"
        run "cp $TMPDIR/$filteredBam.md5 ." "$ECHO"
        log "done\n"
        log "Copying filtered bam file to mapping dir..." $step
        run "cp $TMPDIR/$filteredBam ." "$ECHO"
        log "done\n"
    #else
    #    log "Error producing filtered bam file" "ERROR" >&2
    #    exit -1
    fi
    endTime=$(date +%s)
    printHeader "Conversion step completed in $(echo "($endTime-$startTime)/60" | bc -l | xargs printf "%.2f\n") min"
else
    printHeader "Bam file is present...skipping conversion step"
fi

pipelineEnd=$(date +%s)

log "\n"
printHeader "Blueprint pipeline (only bam file) for $sample completed in $(echo "($pipelineEnd-$pipelineStart)/60" | bc -l | xargs printf "%.2f\n") min "
exit 0
