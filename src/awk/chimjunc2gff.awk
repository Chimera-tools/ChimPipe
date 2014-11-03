#!/usr/bin/env awk

# *****************************************************************************

#       chimjunc2gff.awk

#       This file is part of the ChimPipe pipeline 

#       Copyright (c) 2014 Bernardo Rodríguez-Martín 
#                                          Emilio Palumbo 
#                                          Sarah djebali 

#       Computational Biology of RNA Processing group
#       Department of Bioinformatics and Genomics
#       Centre for Genomic Regulation (CRG)
                                           
#       Github repository - https://github.com/Chimera-tools/ChimPipe

#       Documentation - https://chimpipe.readthedocs.org/

#       Contact - chimpipe.pipeline@gmail.com

#       Licenced under the GNU General Public License 3.0 license.
#******************************************************************************

# Description
##############
# Convert a chimeric junctions into gff format. From each chimeric junction produces two GFF features, one for each side of the junction. 
# Deals with stranded and unstranded data

## Input: Chimeric jucntions
# chr7_128502985_+:chr1_231556974_- 1 4 128502926 231556959
# 10777 (5 fields)   *** juncid, nbstag, nbtot, beg, end (from staggered_to_junction.awk)

## Output: GFF
# chr7  . exon 128502926 128502985 . + . junc: chr7_128502985_+:chr1_231556974_-
# chr11 . exon 231556959 231556974 . - . junc: chr7_128502985_+:chr1_231556974_-

# Usage example:
# awk -v stranded=$stranded -f chimjunc2gff.awk

# Where $stranded can be:
# 	 	
# 		0: unstranded
#		1: stranded

{
    split($1,a,":"); 
    split(a[1],a1,"_"); 
    split(a[2],a2,"_"); 
    
    if (stranded==1) # stranded data 
    {	
    	if (a1[2]<$4)	# part 1
    	{
			print a1[1], ".", "exon", a1[2], $4, ".", a1[3], ".", "junc:", $1; 
    	}
    	else
    	{
			print a1[1], ".", "exon", $4, a1[2], ".", a1[3], ".", "junc:", $1;     
    	}
    
    	if (a2[2]<$5)	# part 2
    	{
			print a2[1], ".", "exon", a2[2], $5, ".", a2[3], ".", "junc:", $1; 
    	}
    	else
    	{
			print a2[1], ".", "exon", $5, a2[2], ".", a2[3], ".", "junc:", $1;     
    	}
    }
    
    else  # unstranded data 
    {	
    	if (a1[2]<$4)	# part 1
    	{
			print a1[1], ".", "exon", a1[2], $4, ".", ".", ".", "junc:", $1; 
    	}
    	else
    	{
			print a1[1], ".", "exon", $4, a1[2], ".", ".", ".", "junc:", $1;     
    	}
    
    	if (a2[2]<$5)	# part 2
    	{
			print a2[1], ".", "exon", a2[2], $5, ".", ".", ".", "junc:", $1; 
    	}
    	else
    	{
			print a2[1], ".", "exon", $5, a2[2], ".", ".", ".", "junc:", $1;     
    	}
    }
}    
    


