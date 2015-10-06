
#!/usr/bin/env awk

# *****************************************************************************
	
#	
	
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

### output


# Usage example:
################

BEGIN{
	while (getline < fileRef > 0)
	{
		spliceJunc=$1; 
		geneId=$19; 
		geneName=$21; 
		readIds=$NF; 
		split(readIds,readIdsList,","); 
		
		for (readId in readIdsList)
		{
			mateId=readIdsList[readId];
			
			if (geneIds[mateId] != "unannotated")
			{
				spliceJunctions[mateId]=spliceJunc; 
				geneIds[mateId]=geneId; 
				geneNames[mateId]=geneName;
			}
		}
	}
}
{
	readId=$1;
	gnListA=$5; 
	
	
	split(readId,readIdList,"/"); 
	pairId=readIdList[1]; 
	contiguousMate=readIdList[2];
	
	if (contiguousMate =="1")
	{
		mateId=pairId"/2";
	}
	else
	{	
		mateId=pairId"/1";
	};
		
	# Set DISCORDANT pair as default
	type="DISCORDANT";
		
	
	if (spliceJunctions[mateId]!="")
	{
		gnListB=geneIds[mateId]; 
		nbA=split(gnListA, a, ","); 
		nbB=split(gnListB,b, ","); 
		
		for (i=1;i<=nbA;i++)
		{
			for (j=1;j<=nbB;j++)
			{
				if (a[i]==b[j])
				{
					type="CONCORDANT";
				}
			}
		}
	}
	else
	{
		type="UNPAIRED"
	} 
	
	if (type == "DISCORDANT")
	{
		if (contiguousMate =="1")
		{
			mapCoords1=$3;
			percOverlap1=$4;
			geneId1=$5;
			geneName1=$6;
		
			mapCoords2=spliceJunctions[mateId];
			percOverlap2="na";
			geneId2=geneIds[mateId];
			geneName2=geneNames[mateId];	
		}
		else
		{
			mapCoords1=spliceJunctions[mateId];
			percOverlap1="na";
			geneId1=geneIds[mateId];
			geneName1=geneNames[mateId];
			
			mapCoords2=$3;
			percOverlap2=$4;
			geneId2=$5;
			geneName2=$6;
		}; 
		
		row=pairId"\t"type"\t"mapCoords1"\t"mapCoords2"\t"percOverlap1"\t"percOverlap2"\t"geneId1"\t"geneId2"\t"geneName1"\t"geneName2;
		
		print row;	
	}
}
