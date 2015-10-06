#!/usr/bin/env awk

# *****************************************************************************
	
#	select_overlappingExons_spliceSide.awk
	
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
# chr7	ChimPipe	spliceJunc1	100859976	100860018	.	-	.	JuncId: "chr7_100859976_-:chr7_100859827_-";	chr7	HAVANA	exon	100859976	100860067	.	-	.	gene_name "PLOD3"; gene_id "ENSG00000106397.7"; gene_type "protein_coding"; transcript_id "ENST00000223127.3";	43 100 0
# chr7	ChimPipe	spliceJunc1	100859976	100860018	.	-	.	JuncId: "chr7_100859976_-:chr7_100859827_-";	chr7	HAVANA	exon	100859976	100860123	.	-	.	gene_name "PLOD3"; gene_id "ENSG00000106397.7"; gene_type "protein_coding"; transcript_id "ENST00000489927.1";	43 100 0
# chr17	ChimPipe	spliceJunc1	42766918	42766967	.	-	.	JuncId: "chr17_42766918_-:chr17_42761327_-";	chr17	HAVANA	exon	42766918	42767093	.	-	.	gene_name "CCDC43"; gene_id "ENSG00000180329.9"; gene_type "protein_coding"; transcript_id "ENST00000588687.1";	50 100 0
# chr17	ChimPipe	spliceJunc1	42766918	42766967	.	-	.	JuncId: "chr17_42766918_-:chr17_42761327_-";	chr17	HAVANA	exon	42766918	42767113	.	-	.	gene_name "CCDC43"; gene_id "ENSG00000180329.9"; gene_type "protein_coding"; transcript_id "ENST00000592333.1";	50 100 0

## output
# chr7_100859976_-:chr7_100859827_- spliceJunc1 100 0 ENSG00000106397.7:chr7_100859976_100859827_- ENSG00000106397.7 PLOD3 protein_coding
# chr17_42766918_-:chr17_42761327_- spliceJunc2 100 0 ENSG00000180329.9:chr17_42766918_42761327_- ENSG00000180329.9 CCDC43 protein_coding


# Usage example:
################
# awk -f select_overlappingExons_spliceSide.awk file.txt

# Make empty vectors to save information about each splice junction block: 
BEGIN{
	split("", gnIds); # ids of the different genes the junction block overlaps
	split("", gnNames);
	split("", gnTypes); 
	split("", percOvelaps); 
	split("", distExonSS); 
}

{
	# Initialize variables to default values
	overlap="0";
	dist="na";
	gnId="unannotated"; 
	gnName="na"; 
	gnType="na"; 
	 
	addGn="1";
	
	# Check which block of the splice junction is and make variable with its coordinate and the block identifier
	blockSide=$3;	
	
	split($10,a,"\""); 
	juncCoord=a[2];
	
	
	juncBlockId=blockSide":"juncCoord
	juncCoords[juncBlockId]=juncCoord
	
	
	# Exon coordinates, percentage of overlap between junction block and a given exon and distance between the junction coordinates and the exon splice site 
	if ($NF != "-1")
	{
		overlap=$(NF-1);
		dist=$NF;
	}
	
	# Parse exon atributes list and select if available its gene id, name and type	
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
		else if ($i=="gene_type")
		{
			split($(i+1),type,"\""); 
			gnType=type[2];
		}
	};
		
	n1=split(gnIds[juncBlockId],ids,","); 
	
	# A) Better choice -> Junction block do not overlapping any exon OR if empty list (no previous overlapping exon for this junction side)  
	# OR if the distance to the exon boundary is lower than the previous one OR the percentage of overlap 
	# is higher than the previous one
	if ((gnId=="unAnnotated")||(n1 == "0")||((dist <= exonBoundDist[juncBlockId])&&(overlap >= percOvelaps[juncBlockId])&&((dist != exonBoundDist[juncBlockId])||(overlap != percOvelaps[juncBlockId]))))
	{
		gnIds[juncBlockId]=gnId; 
		gnNames[juncBlockId]=gnName;
		gnTypes[juncBlockId]=gnType;  
		exonBoundDist[juncBlockId]=dist;
		percOvelaps[juncBlockId]=overlap;
	
		addGn="0";
	}
	# B) At least equally good choice <- Junction block overlapping an exon with a percentage of overlap 
	# AND a distance to the exon boundary equal to the previous one
	else if ((dist == exonBoundDist[juncBlockId])&&(overlap == percOvelaps[juncBlockId]))
	{	
		n2=split(gnId, gnIdsExon, ",");
		
		# Parse list of genes overlaping the junction block
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
	# C) Worse choice, percentage of overlap lower than the previous one and/or 
	# the distance to the exon boundary is higher 
	else 
	{
		addGn="0";	
	}
	
	
	# Add the gnId, name and biotype to the list
	if (addGn=="1")
	{
		gnIds[juncBlockId]=gnId","gnIds[juncBlockId]; 
		gnNames[juncBlockId]=gnName","gnNames[juncBlockId];
		gnTypes[juncBlockId]=gnType","gnTypes[juncBlockId];  
	}
}

END{
	for (block in juncCoords)
	{
		split(block, a, ":");
		blockSide=a[1];

		print juncCoords[block], blockSide, percOvelaps[block], exonBoundDist[block], gnIds[block], gnNames[block], gnTypes[block];
	}
}

