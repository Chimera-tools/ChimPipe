#!/usr/bin/env awk

# *****************************************************************************
	
#	bed2gff.awk
	
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
# Convert BED record with only one "blockCount" into gff format.

## Input: BED 
#chr1    12117   12170   SRR201779.1488862_PATHBIO-SOLEXA2_30TUEAAXX:3:22:1724:187_length=53#0/1 1       -       12117   12170   255,0,0 1       53      0
#chr1    12154   12207   SRR201779.673684_PATHBIO-SOLEXA2_30TUEAAXX:3:10:1411:1581_length=53#0/2 8       -       12154   12207   255,0,0 1       53      0

## Output: GFF
# chr1 ChimPipe alblock 12118 12170 1 - . name: SRR201779.1488862_PATHBIO-SOLEXA2_30TUEAAXX:3:22:1724:187_length=53#0/1
# chr1 ChimPipe alblock 12155 12207 8 - . name: SRR201779.673684_PATHBIO-SOLEXA2_30TUEAAXX:3:10:1411:1581_length=53#0/2


BEGIN{OFS="\t"}
{
	block=$1"\tChimPipe\talBlock\t"($2+1)"\t"$3"\t"$5"\t"$6"\t.\tReadName: \""$4"\"\;";
	print block;
}


