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
# awk -v OFS="\t" <INPUT_ARGS> -f SAMfilter.awk


# Where INPUT_ARGS can be:
# -v higherThan="1". 1: print multimapped and uniquely mapped reads with a number of multimappings and mismatches, respectivelly, higher than a threshold; 0: print multimapped and uniquely mapped reads with a number of multimappings and mismatches, respectivelly, equal or lower than a threshold

# -v unmapped="1"  Print unmapped reads;  
# -v unique="1" Print uniquelly mapped reads;
# -v multimapped="1" Print multi-mapped reads;
# -v nbMatches="N" Number of multimappings threshold
# -v nbMism="N" Number of mismatches threshold
# -v NHfield="" Position of the field containing the number of alignments (NH). Mandatory if "unique=1" or "multimapped=1"
# -v NMfield="" Position of the field containing the edit distance (NM).
# -v higherThan="1" Print uniquely mapped reads with a number of mismatches higher than a treshold or multimapped with a number of hits higher than treshold. If not set, do the opposite, comparison operator lower than or equal instead of higher than.

## Examples:

# Extract reads uniquely mapped with 4 or less mismatches:
# awk -v OFS="\t" -v unique="1"  -v nbMism="4" -v NHfield="13" -v NMfield="14" -f SAMfilter.awk 

# Extract reads Unmapped, multimapped and uniquely mapped with more than 4 mismatches: 
# awk -v OFS="\t" -v higherThan="1" -v unmapped="1" -v multimapped="1" -v unique="1"  -v nbMism="4" -v NHfield="13" -v NMfield="14" -f SAMfilter.awk 
 
 
BEGIN{
	if	(((unique=="1")||(multimapped=="1"))&&(NHfield==""))
	{
		print "[ERROR] Please specify the position of the NH field (nb. of hits). I.e: awk -v NHfield=14 ..." > "/dev/stderr";
    	exit 1
	}
	else
	{
		if ((unique=="1")&&(nbMism!="")&&(NMfield==""))
		{
			print "[ERROR] Please specify the position of the NM field (edit distance). I.e: awk -v NMfield=15 ..." > "/dev/stderr";
    		exit 1
		} 
	}	
}
 
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
	else
	{
		if (higherThan=="1") # Flag to set the nbMism and nbMatches comparison operators to > enabled
		{
			if ((unique=="1")&&(NH[3]=="1")&&((nbMism=="")||(NM[3]>nbMism))) # Print all uniquelly mapped reads or unique with more than X mismatches if "nbMism" input variable set
			{
					print $0;
			}
			else if ((multimapped=="1")&&(NH[3]>"1")&&((nbMatches=="")||(NH[3]>nbMatches))) # Print all multimapped reads or multimappings with more than X hits if "nbMatches" input variable set
			{
					print $0;
			}
		}
		else  # By default the comparison operator is <=
		{
			if ((unique=="1")&&(NH[3]=="1")&&((nbMism=="")||(NM[3]<=nbMism))) # Print all uniquelly mapped reads or unique with X or less mismatches if "nbMism" input variable set
			{
				print $0;
			}
			else if ((multimapped=="1")&&(NH[3]>"1")&&((nbMatches=="")||(NH[3]<=nbMatches))) # Print all multimapped reads or multimappings with X or less hits if "nbMatches" input variable set
			{
				print $0;
			}
		}	
	}
}

