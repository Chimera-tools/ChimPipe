#!/bin/bash

if [[ $# != 2 ]]; then
    echo "Usage: infer_exp.sh ANNOTATION_FILE BAM_FILE"
    exit 1
fi
# binaries
gtfToGenePred=/users/rg/epalumbo/bin/gtfToGenePred
genePredToBed12=/users/rg/epalumbo/bin/genePredToBed12.awk

# I/O
anno=$1
bam=$2
genePred=`basename $anno .gtf`.genePred
bed12=`basename $anno .gtf`.bed

echo "Creating the reference gene model in bed format" >&2
if [ ! -e $genePred ]; then
    $gtfToGenePred $anno -allErrors $genePred 2> $genePred.err
fi
if [ ! -e $bed12 ]; then
    cat $genePred | $genePredToBed12 > $bed12
fi

echo "Inferring experiment" >&2
. /software/rg/el6.3/virtualenvs/python2.7.3/bin/activate
infer_experiment.py -i $bam -r $bed12
echo "Removing temporary files" >&2
rm $genePred $genePred.err $bed12
echo "DONE" >&2
deactivate
