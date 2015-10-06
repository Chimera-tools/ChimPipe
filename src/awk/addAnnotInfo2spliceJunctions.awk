#!/usr/bin/env awk

# *****************************************************************************
	
#	addAnnotInfo2spliceJunctions.awk
	
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
# Takes as input a file containing the annotated splice junction sides and another one with the supporting reads among other information and produces a combined output file. 

## input 1 (fileRef)
# chr7_100859976_-:chr7_100859827_- spliceJunc1 100 0 ENSG00000106397.7 PLOD3 protein_coding

## input 2 

## output
# chr7_100859976_-:chr7_100859827_-	10	11	11	0	100860018	100859780	GT	AG	1	1	149	100	100	0	0ENSG00000106397.7	ENSG00000106397.7	PLOD3	PLOD3	protein_coding	protein_coding	SRR201779.3750664_PATHBIO-SOLEXA2_30TUEAAXX:3:58:1707:1197_length=53#0/1,SRR201779.5891847_PATHBIO-SOLEXA2_30TUEAAXX:3:91:769:114_length=53#0/1,SRR201779.3883456_PATHBIO-SOLEXA2_30TUEAAXX:3:61:823:1701_length=53#0/2,SRR201779.1222929_PATHBIO-SOLEXA2_30TUEAAXX:3:18:21:810_length=53#0/1,SRR201779.4428441_PATHBIO-SOLEXA2_30TUEAAXX:3:69:484:1727_length=53#0/1,SRR201779.3750664_PATHBIO-SOLEXA2_30TUEAAXX:3:58:1707:1197_length=53#0/2,SRR201779.5219014_PATHBIO-SOLEXA2_30TUEAAXX:3:81:330:424_length=53#0/1,SRR201779.6481993_PATHBIO-SOLEXA2_30TUEAAXX:3:99:1639:590_length=53#0/1,SRR201779.1851231_PATHBIO-SOLEXA2_30TUEAAXX:3:28:1220:256_length=53#0/2,SRR201779.337899_PATHBIO-SOLEXA2_30TUEAAXX:3:5:1084:1720_length=53#0/2,SRR201779.3760532_PATHBIO-SOLEXA2_30TUEAAXX:3:58:507:1518_length=53#0/2,


# Usage example:
################
# awk -v juncAnnotated=file1.txt -f addAnnotInfo2spliceJunctions.awk file2.txt

BEGIN{
	## Make a hash for each junction side containing its annotation info
	while (getline < juncAnnotated > 0)
	{
		juncCoord=$1; 
		block=$2; 
		blockId=juncCoord":"block; 
		overlaps[blockId]=$3; 
		distExonBoundaries[blockId]=$4; 
		gnIds[blockId]=$5;
		gnNames[blockId]=$6; 
		gnTypes[blockId]=$7;
	}
}

{	
	juncCoord=$1; 
	nbStag=$2; 
	nbTotal=$3; 
	nbUnique=$4; 
	nbMulti=$5; 
	donorSS=$6; 
	acceptorSS=$7;
	beg=$8; 
	end=$9; 
	sameChrStr=$10;
	okGxOrder=$11;
	dist=$12; 
	readIds=$13; 

	# Gather the annotation info for the splice junction donor side
	overlapA=overlaps[juncCoord":spliceJunc1"]; 
	distExonBoundaryA=distExonBoundaries[juncCoord":spliceJunc1"]; 
	gnIdA=gnIds[juncCoord":spliceJunc1"];
	gnNameA=gnNames[juncCoord":spliceJunc1"]; 
	gnTypeA=gnTypes[juncCoord":spliceJunc1"]; 
	
	# Gather the annotation info for the splice junction acceptor side
	overlapB=overlaps[juncCoord":spliceJunc2"]; 
	distExonBoundaryB=distExonBoundaries[juncCoord":spliceJunc2"];  
	gnIdB=gnIds[juncCoord":spliceJunc2"]; 
	gnNameB=gnNames[juncCoord":spliceJunc2"]; 
	gnTypeB=gnTypes[juncCoord":spliceJunc2"];       
	   	 
	# Compute the percentage of staggered and multimapped reads supporting a junctions
	percStag=nbStag/nbTotal*100;
	percMulti=nbMulti/nbTotal*100;
	
	# Print splice junction plus associated info as output   	 
	
	   	 row=juncCoord"\t"nbTotal"\t"nbStag"\t"percStag"\t"nbMulti"\t"percMulti"\t"overlapA"\t"overlapB"\t"distExonBoundaryA"\t"distExonBoundaryB"\t"donorSS"\t"acceptorSS"\t"beg"\t"end"\t"sameChrStr"\t"okGxOrder"\t"dist"\t"gnIdA"\t" gnIdB"\t"gnNameA"\t"gnNameB"\t"gnTypeA"\t"gnTypeB"\t"readIds;
	print row;  
}
