#!/usr/bin/env awk

# *****************************************************************************
	
#	bedPE2gff.awk
	
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
###############
# Convert BED paired-end into gff format.

## Input: BEDPE
# chr1 752991 753018 chr1 754245 754250 HWI-ST985:73:C08BWACXX:8:2208:2017:40383/1  1000 + +

## Output: GFF
# chr1 ChimPipe  752992 753018 . + . name: HWI-ST985:73:C08BWACXX:8:2208:2017:40383/1
# chr1 ChimPipe alblock 754246 754250 . + . name: HWI-ST985:73:C08BWACXX:8:2208:2017:40383/1

BEGIN{OFS="\t"}
{
	block1=$1"\tChimPipe\talBlock1\t"($2+1)"\t"$3"\t"$8"\t"$9"\t.\tReadName: \""$7"\"\;";
	block2=$4"\tChimPipe\talBlock2\t"($5+1)"\t"$6"\t"$8"\t"$10"\t.\tReadName: \""$7"\"\;";
	
	print block1;
	print block2;
}
