#!/usr/bin/env awk

# *****************************************************************************
	
#	makeStaggeredSplitMappings.awk
	
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

# Takes as input a file 

# input
# chr1	ChimPipe	alBlock1	55622718	55622724	1	+	.	ReadName: "SRR064438.3006298_HWI-EAS418:5:46:617:1074_length=50#0/1";
# chr1	ChimPipe	alBlock2	55622929	55622971	1	+	.	ReadName: "SRR064438.3006298_HWI-EAS418:5:46:617:1074_length=50#0/1";
# chr2	ChimPipe	alBlock1	110342896	110342898	3	+	.	ReadName: "SRR064287.3353152_HWI-EAS418:7:49:1321:1518_length=50#0/2";
# chr2	ChimPipe	alBlock2	110343298	110343344	3	+	.	ReadName: "SRR064287.3353152_HWI-EAS418:7:49:1321:1518_length=50#0/2";

# output

# Usage example:
################
# awk  -f makeStaggeredSplitMappings.awk splitMappings.gff


{
	chrA=$1; 
	begA=$4; 
	endA=$5; 
	strandA=$7; 
	nbHits=$6; 
	
	getline; 
	chrB=$1; 
	begB=$4; 
	endB=$5; 
	strandB=$7; 
	
	
	staggered=chrA"_"begA"_"endA"_"strandA":"chrB"_"begB"_"endB"_"strandB; 
	
	split($10,a,"\""); 
	totalNbSupportReads[staggered]++; 
	supportReads[staggered]=a[2]","supportReads[staggered]; 
	
	if (nbHits=="1"){
		nbUniqueSupportReads[staggered]++;
	}
	else
	{
		nbMultiSupportReads[staggered]++;
	}
}
END{
	for (staggered in totalNbSupportReads)
	{
		totalNb=totalNbSupportReads[staggered]; 
		allSupportReads=supportReads[staggered]; 
		nbUnique=(nbUniqueSupportReads[staggered]!="" ? nbUniqueSupportReads[staggered] : "0"); 
		nbMulti=(nbMultiSupportReads[staggered]!="" ? nbMultiSupportReads[staggered] : "0");  
		
		print staggered, totalNb, nbUnique, nbMulti, allSupportReads; 
	}
}
