#!/bin/bash

<<authors
*****************************************************************************
	
	detect.fq.qual.sh
	
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

function usage
{
cat <<help
	Usage:  detect.fq.qual.sh file.fastq[.gz]  
    
       Takes as input fastq file that can be gzipped or not and output Offset followed by the quality offset of the reads for illumina reads
    exit 0
help
}

# will exit if there is an error or in a pipe
set -e -o pipefail

# GETTING INPUT ARGUMENTS 
#########################
file=$1

# SETTING VARIABLES AND INPUT FILES
###################################
if [[ ! -e $file ]]; then printf "\n\tERROR: Please specify a valid fastq file\n\n" >&2; usage; exit -1; fi


# START
########
if [[ $file ]]; then
    command="cat"
    if [[ $file =~ .*\.gz ]];then
        command="zcat"
    fi
    command="$command $file | "
fi

command="${command}awk 'BEGIN{for(i=1;i<=256;i++){ord[sprintf(\"%c\",i)]=i;}}NR%4==0{split(\$0,a,\"\");for(i=1;i<=length(a);i++){if(ord[a[i]]<59){print \"Offset 33\";exit 0}if(ord[a[i]]>74){print \"Offset 64\";exit 0}}}'"

eval $command

