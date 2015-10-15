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
# Mate 1 file
# SRR201779.4894699_PATHBIO-SOLEXA2_30TUEAAXX:3:76:10:1139_length=53#0/1  chr1_565366_565418_+    ENSG00000225630.1
# ..

# Mate2 file
# SRR201779.4894699_PATHBIO-SOLEXA2_30TUEAAXX:3:76:10:1139_length=53#0/2  chr1_150280531_150280583_-      ENSG00000187145.10
# ..

### output
# SRR201779.4894699_PATHBIO-SOLEXA2_30TUEAAXX:3:76:10:1139_length=53#0    DISCORDANT      chr1_565366_565418_+    chr1_150280531_150280583_-      ENSG00000225630.1       ENSG00000187145.10
# SRR201779.5549320_PATHBIO-SOLEXA2_30TUEAAXX:3:86:252:328_length=53#0    CONCORDANT      chr1_17483_17535_-      chr1_17321_17373_+      ENSG00000227232.4       ENSG00000227232.4
# SRR201779.4899422_PATHBIO-SOLEXA2_30TUEAAXX:3:76:831:388_length=53#0    UNANNOTATED     chr1_17606_17658_-      chr1_17444_17496_+      ENSG00000227232.4       unannotated
# SRR201779.2689262_PATHBIO-SOLEXA2_30TUEAAXX:3:41:1598:622_length=53#0/1 UNPAIRED        chr1_139561_139613_+    unannotated


# Usage example:
################

BEGIN{
	## print header
	header="pairId\ttype\tmapCoordsA\tmapCoordsB\tgnIdA\tgnIdB";	
	print header
	
	while (getline < fileRef > 0)
	{
		readId=$1;
		split(readId,readIdList,"/"); 
		pairId=readIdList[1]; 
	
		mapCoords2=$2;
		geneId2=$3;
		
		mate2[pairId]=mapCoords2"|"geneId2"|1";
	}
}
{
	#print "classifying..";
	readId=$1;
	split(readId,readIdList,"/"); 
	pairId=readIdList[1]; 
	
	mapCoords1=$2;		
	geneId1=$3; 
	
	
	# Set DISCORDANT pair as default
	type="DISCORDANT";
		
	# A) Both mates mapping
	if ( mate2[pairId] != "")
	{
		split(mate2[pairId],mate2List,"|");
		
		mapCoords2=mate2List[1];
		geneId2=mate2List[2];
		
		mate2[pairId]=geneId2"|"mapCoords2"|0";
		
		# A.A) Both mates mapping and at least one of them do not overlap any annotated gene 
		if ((geneId1 == "unannotated") || (geneId2 == "unannotated"))
		{
			type="UNANNOTATED";
		}
		# C) Both mates mapping and both overlap an annotated gene
		else
		{
		    nb1=split(geneId1, gnList1, ","); 
		    nb2=split(geneId2, gnList2, ",");
			
		    # Do all the possible gene pairs comparisons between the 2 list of genes overlapping each mate. 
		    # If both mates overlap the same gene, set type to concordant. 
		    for (i=1;i<=nb1;i++)
		    {
				for (j=1;j<=nb2;j++)
				{
			    	if (gnList1[i]==gnList2[j])
			    	{
						type="CONCORDANT";
			    	}
				}
		   	} 
		}
		
		row=pairId"\t"type"\t"mapCoords1"\t"mapCoords2"\t"geneId1"\t"geneId2;
	}
	# B) Mate two do not mapping
	else
	{
		type="UNPAIRED";
		row=readId"\t"type"\t"mapCoords1"\t"geneId1;
	}
		
	print row;
}
END{
	for (pairId in mate2) 
	{	
		split(mate2[pairId],mate2List,"|");
		
		boolean=mate2List[3];
		
		if (boolean == "1" )
		{
			mapCoords2=mate2List[1];
			geneId2=mate2List[2];
		
			type="UNPAIRED";
			readId=pairId"/2";
			
			row=readId"\t"type"\t"mapCoords2"\t"geneId2;
			
			print row;			
		}
	}
		
} 

