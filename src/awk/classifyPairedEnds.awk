#!/usr/bin/env awk

# *****************************************************************************
	
#	classifyPairedEnds.awk
	
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
# ........

### input
# SRR201779.3622811_PATHBIO-SOLEXA2_30TUEAAXX:3:56:1614:1647_length=53#0/1 100 ENSG00000184012.7:chr21_42842575_42842670 ENSG00000184012.7 TMPRSS2 protein_coding
# SRR201779.243167_PATHBIO-SOLEXA2_30TUEAAXX:3:4:1747:316_length=53#0/1 100 ENSG00000184012.7:chr21_42836480_42838080,ENSG00000184012.7:chr21_42836478_42838080 ENSG00000184012.7 TMPRSS2 protein_coding

### output


# Usage example:
################

BEGIN{
	## print header

	header="pairId\ttype\tmapCoordsA\tmapCoordsB\tgnIdA\tgnIdB\tgnNameA\tgnNameB\toverlapA\toverlapB";	
	print header
}
{
	split($1,readId,"/"); 
	
	pairId=readId[1]; 
	mate=readId[2]; 


	if (mate=="1")
	{
		#print "eii 1", $0;
		mapCoords1[pairId]=$2;		
		percOverlap1[pairId]=$3; 
		geneId1[pairId]=$4; 
		geneName1[pairId]=$5; 
		pairIds[pairId]="1";
	}
	else
	{
		#print "eii 2", $0;
		mapCoords2[pairId]=$2;
		percOverlap2[pairId]=$3; 
		geneId2[pairId]=$4; 
		geneName2[pairId]=$5;
		pairIds[pairId]="1";
	}
}
END{
	for (pairId in pairIds)
	{		
		# Set DISCORDANT pair as default
		type="DISCORDANT";
		
		# A) One of the mates is not mapping 
		if ((geneId1[pairId]=="") || (geneId2[pairId]=="")) 
		{
			if (geneId1[pairId]=="")
			{
				mateId=pairId"/2";
				# print "mate 2", mateId;
			}
			else
			{
				mateId=pairId"/1";
				# print "mate 1", mateId;
			}
			type="UNPAIRED";
		}	
		# B) Both mates mapping and at least one of them do not overlap any annotated gene 
		else if ((geneId1[pairId] == "unannotated") || (geneId2[pairId] == "unannotated"))
		{	
			type="UNANNOTATED";
		}
		# C) Both mates mapping and both overlap an annotated gene
		else
		{
		    # If gene name available use gene name (better, since for gencode there are several different gene ids for some gene within the same annotation)
		    if ((geneName1[pairId]!="na")&&(geneName2[pairId]!="na"))
		    {
		    	nb1=split(geneName1[pairId], gnList1, ","); 
		    	nb2=split(geneName2[pairId], gnList2, ",");
			}
			# If not, use gene id
			else
			{
				nb1=split(geneId1[pairId], gnList1, ","); 
		    	nb2=split(geneId2[pairId], gnList2, ",");
			}
			
		    # Do all the possible gene pairs comparisons between the 2 list of genes overlapping each mate. 
		    # If both mates overlap the same gene, set type to concordant. 
		    for (i=1;i<=nb1;i++)
		    {
				for (j=1;j<=nb2;j++)
				{
			    	#print gnList1[i],  gnList2[j]; 
			
			    	if (gnList1[i]==gnList2[j])
			    	{
						type="CONCORDANT";
			    	}
				}
		   	} 
		}

		# Report classified read pairs 
		if (type=="UNPAIRED")
		{
			if (mapCoords1[pairId] != "")
			{
				row=mateId"\t"type"\t"mapCoords1[pairId]"\t"percOverlap1[pairId]"\t"geneId1[pairId]"\t"geneName1[pairId];
			}
			else
			{
				row=mateId"\t"type"\t"mapCoords2[pairId]"\t"percOverlap2[pairId]"\t"geneId2[pairId]"\t"geneName2[pairId];
			}
			
		}
		else
		{
			row=pairId"\t"type"\t"mapCoords1[pairId]"\t"mapCoords2[pairId]"\t"percOverlap1[pairId]"\t"percOverlap2[pairId]"\t"geneId1[pairId]"\t"geneId2[pairId]"\t"geneName1[pairId]"\t"geneName2[pairId];
		}
		
		print row;
	}
} 

