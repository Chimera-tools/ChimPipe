#!/usr/bin/env awk

# *****************************************************************************
	
#	select_overlappingExons_contiguousMapping.awk
	
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

# Takes as input a file produced with "compute_juncPercOverlap2Exon_juncDist2ExonSS.awk" For each splice junction side select those overlaping exons which maximize the percentage of overlap and minimize the distance to the exon boundary

## input
#chr1	ChimPipe	alBlock	11856	11908	8	+	.	ReadName: "SRR201779.3905018_PATHBIO-SOLEXA2_30TUEAAXX:3:61:1631:419_length=53#0/1";	chr1	HAVANA	exon	11869	12227	.	+	.	gene_name "DDX11L1"; gene_id "ENSG00000223972.4"; gene_type "pseudogene"; transcript_id "ENST00000456328.2";	40 75.4717

## output


# Usage example:
################
# awk -f select_overlappingExons_contiguousMapping.awk file.txt

# Make empty vectors to save information about each contiguous mapping block: 
BEGIN{
	split("", gnIds); # ids of the different genes mapping block overlaps
	split("", gnNames);
	split("", percOvelaps);  
}

{
	# Initialize variables to default values
	overlap="0";
	gnId="unannotated"; 
	gnName="na"; 
	 
	addGn="1";
	
	# Set the mapping block id as the read name
	split($10,a,"\""); 
	readId=a[2];
	
	# Read mapping coordinates
	mapCoord=$1"_"$4"_"$5"_"$7;
	
	# Percentage of overlap between mapping block and a given exon 
	if ($NF != "-1")
	{
		overlap=$NF;
	}
	
	# Parse exon atributes list and select if available its gene id and name 	
	for(i=19; i<=(NF); i++)
	{
		if ($i=="gene_id")
		{
			split($(i+1),id,"\""); 
			gnId=id[2]; 
		}
		else if ($i=="gene_name") 
		{
			split($(i+1),name,"\""); 
			gnName=name[2];
		}		
	};
	
	n1=split(gnIds[readId],ids,","); 
	
	
	# A) Better choice -> mapping block do not overlapping any exon OR if empty list (no previous overlapping exon for this mapping block)  
	# OR the percentage of overlap is higher than the previous one
	if ((gnId=="unannotated") || (n1 == "0") || (overlap > percOvelaps[readId]))
	{
		mapCoords[readId]=mapCoord;
		percOvelaps[readId]=overlap;
		gnIds[readId]=gnId;   
		gnNames[readId]=gnName;
		
	
		addGn="0";
	}
	# B) At least equally good choice <- Mapping block overlapping an exon with a percentage of overlap equal to the previous one 
	else if (overlap == percOvelaps[readId])
	{	
		n2=split(gnId, gnIdsExon, ",");
		
		# Parse list of genes overlaping the read alignment
		for (i=1; i<=n1; i++)
		{	
			for (j=1; j<=n2; j++)
			{								
				# Do not add the gene id if the current overlapping exon belongs to a gene already in the list 
				if (ids[i]==gnIdsExon[j])
				{
					addGn="0";
				}	
			}	
		}
	}
	# C) Worse choice, percentage of overlap lower than the previous one. 
	else 
	{
		addGn="0";	
	}
	
	# Add the gnId, name and biotype to the list
	if (addGn=="1")
	{
		gnIds[readId]=gnId","gnIds[readId];   
		gnNames[readId]=gnName","gnNames[readId];
	}
}

END{
	for (readId in gnIds)
	{
		print readId, mapCoords[readId], percOvelaps[readId], gnIds[readId], gnNames[readId];
	}
}

