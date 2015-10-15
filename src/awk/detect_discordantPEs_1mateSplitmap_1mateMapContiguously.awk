
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
## Chimsplice normal junctions
# chr3_182637125_+:chr3_182637126_+       NORMAL  1       1       100     0       0       100     100     2298    1315    GT      AG      182637098       182637150       1       1       1   ENSG00000058063.11       ENSG00000058063.11      ATP11B  ATP11B  protein_coding  protein_coding  SRR201779.3607803_PATHBIO-SOLEXA2_30TUEAAXX:3:56:1701:1083_length=53#0/1,
# chr4_71629386_+:chr4_71630192_+ NORMAL  2       2       100     0       0       100     100     0       0       GT      AG      71629354        71630235        1       1       806     ENSG00000018189.8    ENSG00000018189.8       RUFY3   RUFY3   protein_coding  protein_coding  SRR201779.54315_PATHBIO-SOLEXA2_30TUEAAXX:3:1:1271:382_length=53#0/1,SRR201779.1335758_PATHBIO-SOLEXA2_30TUEAAXX:3:20:515:1865_length=53#0/1,
# chr22_24563281_+:chr22_24564400_+       NORMAL  2       2       100     0       0       100     68.0851 0       15      GT      AG      24563252        24564446        1       1       1119ENSG00000099991.12       ENSG00000099991.12      CABIN1  CABIN1  protein_coding  protein_coding  SRR201779.3315121_PATHBIO-SOLEXA2_30TUEAAXX:3:51:533:1902_length=53#0/2,SRR201779.5842695_PATHBIO-SOLEXA2_30TUEAAXX:3:90:1461:1795_length=53#0/1,

## Unpaired
# SRR201779.2689262_PATHBIO-SOLEXA2_30TUEAAXX:3:41:1598:622_length=53#0/1 UNPAIRED        chr1_139561_139613_+    unannotated
# SRR201779.5987634_PATHBIO-SOLEXA2_30TUEAAXX:3:92:1029:828_length=53#0/1 UNPAIRED        chr1_169214_169266_+    ENSG00000241860.2
# SRR201779.5090138_PATHBIO-SOLEXA2_30TUEAAXX:3:79:1575:982_length=53#0/1 UNPAIRED        chr1_565178_565228_+    ENSG00000225630.1
# SRR201779.2460298_PATHBIO-SOLEXA2_30TUEAAXX:3:37:508:1042_length=53#0/1 UNPAIRED        chr1_565294_565346_+    ENSG00000225630.1

### output


# Usage example:
################

BEGIN{
	while (getline < fileRef > 0)
	{
		spliceJunc=$1; 
		geneId=$19; 
		readIds=$NF; 
		split(readIds,readIdsList,","); 
		
		for (readId in readIdsList)
		{
			mateId=readIdsList[readId];
			
			if (geneIds[mateId] != "unannotated")
			{
				spliceJunctions[mateId]=spliceJunc; 
				geneIds[mateId]=geneId; 
			}
		}
	}
}
{
	readId=$1;
	gnListA=$4; 
	
	
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
			geneId1=$4;
		
			mapCoords2=spliceJunctions[mateId];
			geneId2=geneIds[mateId];
		}
		else
		{
			mapCoords1=spliceJunctions[mateId];
			geneId1=geneIds[mateId];
			
			mapCoords2=$3;
			geneId2=$4;
		}; 
		
		row=pairId"\t"type"\t"mapCoords1"\t"mapCoords2"\t"geneId1"\t"geneId2;
		
		print row;	
	}
}
