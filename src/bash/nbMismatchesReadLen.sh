#!/bin/bash

<<authors
*****************************************************************************
	
	nbMismatchesReadLen.sh
	
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
	
	nbMismatchesReadLen.sh -  Compute the maximum number of mismatches for mapping filtering based on the average read-pair length and a given percentage of mismatches
	
	Usage: nbMismatchesReadLen.sh file.fastq1[.gz] file.fastq1[.gz] percMism
	
	Example: nbMismatchesReadLen.sh file.fastq1.gz file.fastq1.gz 4
help
}


# In case the user does not provide any input file
###################################################
if [[ ! -e "$1" ]] || [[ ! -e "$2" ]] || [[ "$3" == "" ]]
then
	printf "\n\tERROR: Please specify a valid input \n\n" >&2; 
	usage; 
	exit -1;
fi

fastq1=$1;
fastq2=$2;
percMism=$3;


# Directories 
#############
## Set root directory
path="`dirname \"$0\"`"              # relative path
rootDir="`( cd \"$path\" && pwd )`"  # absolute path

if [ -z "$rootDir" ] ; 
then
  # error; for some reason, the path is not accessible
  # to the script
  log "Path not accessible to the script\n" "ERROR" 
  exit 1  # fail
fi

## Set bash directory 
bashDir=$rootDir

# Programs and scripts
########################

## Bash
avReadLen=$bashDir/average_readLength.sh

# START
########

# 1. Compute mate 1 and 2 average read lengths
################################################

# Mate 1 average read length
avMate1Len=`bash $avReadLen $fastq1`

# Mate 2 average read length
avMate2Len=`bash $avReadLen $fastq2`

# Average read length for both mates
avMatesLen=`printf "$avMate1Len\t$avMate2Len\n" | awk '{mean = ($1 + $2)/2; print mean;}'`

# 2. Set maximum number of mismatches == X% of the pair length 
################################################################
nbMism=`echo $avMatesLen | awk -v percMism=$percMism '{avMatesLen=$1; nbMism=(avMatesLen*percMism)*2/100; print nbMism;}' | python -c "print int(round(float(raw_input())))"`

echo $nbMism

