#!/usr/bin/env awk

# *****************************************************************************
	
#	select_overlappingGene_contiguousMapping.awk
	
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

# 

## input
# SRR201779.1304693_PATHBIO-SOLEXA2_30TUEAAXX:3:20:525:1065_length=53#0/1 chr1_139634_139686_+    unannotated
# SRR201779.5987634_PATHBIO-SOLEXA2_30TUEAAXX:3:92:1029:828_length=53#0/1 chr1_169214_169266_+    ENSG00000241860.2
# SRR201779.431340_PATHBIO-SOLEXA2_30TUEAAXX:3:7:1362:1943_length=53#0/2  chr1_228510_228562_+    ENSG00000228463.4
# SRR201779.431340_PATHBIO-SOLEXA2_30TUEAAXX:3:7:1362:1943_length=53#0/2  chr1_228510_228562_+    ENSG00000241670.2


## output
# SRR201779.1304693_PATHBIO-SOLEXA2_30TUEAAXX:3:20:525:1065_length=53#0/1 chr1_139634_139686_+    unannotated
# SRR201779.5987634_PATHBIO-SOLEXA2_30TUEAAXX:3:92:1029:828_length=53#0/1 chr1_169214_169266_+    ENSG00000241860.2
# SRR201779.431340_PATHBIO-SOLEXA2_30TUEAAXX:3:7:1362:1943_length=53#0/2  chr1_228510_228562_+    ENSG00000228463.4,ENSG00000241670.2

# Usage example:
################
# awk -f 
BEGIN{
	readId="";
}

{
	# A) First line
	if (NR=="1")
	{
		readId=$1;
		mapCoord=$2;
		gnIds=$3;
	} 
	# B) Not first line and exon overlapping a different read than the one in the previous line
	else if ((NR!="1") && ($1!=readId))
	{
		row=readId"\t"mapCoord"\t"gnIds; 	
		print row;
		
		readId=$1;
		mapCoord=$2;
		gnIds=$3;
	}
	# C) Not first line and exon overlapping the same read as the one in the previous line
	else 
	{
		
		newGnIds=$3;
		
		# Iterate over the two lists of gnIds, add novel gnIds  
		split(newGnIds,newGnIdsList,",");		
		split(gnIds,gnList,",")
		
		for (newId in newGnIdsList)
		{
			addGnId="1";
		
			for (id in gnList)
			{
				if (gnList[id] == newGnIdsList[newId])
				{
					addGnId="0";	
				}	
			}
			
			# Add gene id if overlapping exon belongs to a gene not already in the list
			if (addGnId == "1")
			{	
				gnIds=gnIds","newGnIdsList[newId];
			}
		}	
	}
}

# Last line
END{
	row=readId"\t"mapCoord"\t"gnIds; 	
	print row;
}

