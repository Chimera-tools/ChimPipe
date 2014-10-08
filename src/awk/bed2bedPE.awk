#!/usr/bin/env awk

# *****************************************************************************
	
#	bed2bedPE.awk
	
#	This file is part of the ChimPipe pipeline 

#	Copyright (c) 2014 Bernardo Rodríguez-Martín 
#					   Emilio Palumbo 
#					   Sarah djebali 
	
#	Computational Biology of RNA Processing group
#	Department of Bioinformatics and Genomics
#	Centre for Genomic Regulation (CRG)
					   
#	Github repository - https://github.com/Chimera-tools/ChimPipe
	
#	Documentation - https://chimpipe.readthedocs.org/

#	Contact - chimpipe.pipeline@gmail.com
	
#	Licenced under the GNU General Public License 3.0 license.
#******************************************************************************

# Description
##############
# Convert into BED paired-end format a bed12 file containing split-mapped reads in two blocks. If input variable rev set to different than 0, 
# reverse the blocks when reads split-mapping in -/-

## Input: BED12
#chr22 1000 5000 HWI-ST985:73:C08BWACXX:8:2208:2017:40383/1 960 + 1000 5000 0 2 567,488, 0,3512

## Output: BEDPED
#chr22 1000 1567 chr22 4512 5000 HWI-ST985:73:C08BWACXX:8:2208:2017:40383/1 960 + +

# Usage example:
# awk -v rev="1" -f bed2bedPE.awk

BEGIN{OFS="\t"}
{
	split($11,sz,",");
	split($12,st,",");
	begA=$2+st[1];
	endA=begA+sz[1];
	begB=$2+st[2];
	endB=begB+sz[2];
	if ((rev=="")||(rev==0)||($6=="+")) 
    {
		print $1, begA, endA, $1, begB, endB, $4, $5, $6, $6;
	}
	else
	{
		print $1, begB, endB, $1, begA, endA, $4, $5, $6, $6;
	}
}



