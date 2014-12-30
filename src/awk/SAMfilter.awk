#!/usr/bin/env awk

# *****************************************************************************
	
#	SAMfilter.awk
	
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
# ******************************************************************************

# Description
###############
# Takes as input a SAM file and filters the alignments according to the mapping status: unmapped, multimapped, unique mapped and nb. of mismatches. Filter customizable

##  Input (SAM file)
# SRR201779.3077035_PATHBIO-SOLEXA2_30TUEAAXX:3:47:705:1013_length=53#0	147	chr1	91185	254	30M4I19M	=	91030	-204	AAAGCGGGACTTTCTTAAAAAAAAAAAAAAAAAACTTAGGTTCTTGGATGTTC	.-.)&(&,(&6?93&&.BBBBBCBBBCBBBBCBBCBCBBB?BBBBCCB?BBBC	RG:Z:0	NH:i:1	NM:i:5	XT:A:U	md:Z:19>4-A13G8T4C1                         ** second in pair, read reversed

# SRR201779.1693706_PATHBIO-SOLEXA2_30TUEAAXX:3:25:1514:905_length=53#0	163	chr1	228464	254	53M	=	212482724	212254313	TTGGCTTTCACATCCACCCTGCACCCAAGTGTGTTGTTTACTTCTATCTTCTT	BBBBBBBBBBBA>ABBB=7<BB><@B?9<;&9<BBB;&0&3<A@<B?7?9&93	RG:Z:0	NH:i:1	NM:i:3	XT:A:U	md:Z:8A8AG2A2A4C9T13 ** second in pair, read not reversed

##  Output (Filtered SAM file)

# Usage:
########
# awk -v OFS="\t" -v unmapped=[1|0] -v multimapped=[1|0] -v unique=[1|0] -v higherThan=[1|0] -v nbMism=<NUMBER> -v NHfield=<NUMBER> -v NMfield=<NUMBER>  -f SAMfilter.awk

# Where:
# - unmapped [1|0]. 1: print unmapped reads; 0: discard unmapped reads  
# - multimapped [1|0]. 1: print multi-mapped reads; 0: discard multi-mapped reads
# - unique [1|0]. 1: print uniquely mapped reads with a number of mismatches lower than or equal/higher (controled with "higherThan" parameter) than a treshold  (controled with "nbMism" parameter); 0: discard all the uniquely mapped reads
# - higherThan [1|0]. 1: print reads with a number mismatches higher than a threshold; 0: print reads with a number of mismatches equal or lower than a threshold.  
# - nbMism: number of mismatches threshold
# - NHfield: position of the field containing the number of alignments.
# - NMfield: position of the field containing the edit distance.

## Examples:

# Reads uniquely mapped with 4 or less mismatches:
# awk -v OFS="\t" -v unmapped="0" -v multimapped="0" -v unique="1" -v higherThan="0" -v nbMism="4" -v NHfield="13" -v NMfield="14" -f SAMfilter.awk 

# Reads Unmapped, multimapped and mapped with more than 4 mismatches: 
# awk -v OFS="\t" -v unmapped="1" -v multimapped="1" -v unique="1" -v higherThan="1" -v nbMism="4" -v NHfield="13" -v NMfield="14" -f SAMfilter.awk 
 
 
/^@/{ 	     # Print the header
	print $0;
}

! /^@/{

	split($NHfield,NH,":"); 
	split($NMfield,NM,":"); 
		 
	if ((unmapped=="1")&&(and($2,0x4))) # Print the alignment if unmapped and the unmapped read flag is enabled  
	{
		print $0;
	}	
	else if ((multimapped=="1")&&(NH[3]>1)) # Print the alignment if multimappin and the multimapped read flag is enabled 
	{
		print $0;
	}
	else if ((unique=="1")&&(higherThan=="0")&&(NH[3]=="1")&&(NM[3]<=nbMism)) # Print the alignment if the read maps unique with X or less mismatches
	{
		print $0;
	}
	else if ((unique=="1")&&(higherThan=="1")&&(NH[3]=="1")&&(NM[3]>nbMism))  # Print the alignment if the read maps unique with more than X mismatches
	{
		print $0;
	}
}

