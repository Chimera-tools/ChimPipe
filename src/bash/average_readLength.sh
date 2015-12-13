#!/bin/bash

<<authors
*****************************************************************************
	
	average_readLength.sh
	
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
	
	average_readLength.sh -  Compute the average read length for a fastq file
	
	Usage:  average_readLength.sh file.fastq[.gz]
	
    exit 0
help
}


# GETTING INPUT ARGUMENTS 
#########################
file=$1

# SETTING VARIABLES AND INPUT FILES
###################################
if [[ ! -e $file ]]; then printf "\n\tERROR: Please specify a valid fastq file\n\n" >&2; usage; exit -1; fi

# START
########

if [[ $file =~ .*\.gz ]];
then
	command="zcat";
else
	command="cat";
fi

command="$command $file |"

command="${command}awk 'NR%4==2{readLen=length(\$1); sum=readLen+sum; total++}END{mean = sum/total; print mean;}'"

eval $command


