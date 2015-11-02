#!/usr/bin/env awk

# *****************************************************************************
	
#	extract_kmers_around_splitMappedBlocks.awk
	
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
##############
# Takes as input split-mapped blocks. Select 5-mers coordinates surrounding each of the two mapping blocks. 

BEGIN{

	while (getline < fileRef > 0)
	{
		chrId=$1;
		chrLength=$2;
		
		chrLengths[chrId]=chrLength;
		
		# print chrId, chrLength;
	}
}

{
	#### Get split-mapped block coordinates
	mappingBlocks=$1;
	split(mappingBlocks, a, ":"); 
	
	## Block 1
	split(a[1], a1, "_"); 
	 
	chrA=a1[1]; 
	begA=a1[2]; 
	endA=a1[3]; 
	strA=a1[4]; 
 
 	## Block 2
	split(a[2],a2,"_"); 
	chrB=a2[1]; 
	begB=a2[2]; 
	endB=a2[3]; 
	strB=a2[4]; 
 
 	#### Get 5-mers block coordinates
 	### 1) Block A1
 	
 	## Block beginning
 	# A) 5-mer beginning lower than chromosome beginning coord (lower bound)
 	if ((begA-5) < 1)
 	{
 		begA1="1";
 	}
 	# B) 5-mer beginning within chromosomal coord.
 	else
 	{
 		begA1=(begA-5);
 	}
 	
 	## Block end
 	# A) 5-mer beginning lower than chromosome beginning coord (lower bound)
 	if ((begA-1) < 1)
 	{
 		endA1="1";
 	}
 	# B) 5-mer beginning within chromosomal coord.
 	else
 	{
 		endA1=(begA-1);
 	}
 	
 	## Block coordinates 	
 	blockA1=chrA"\t"strA"\t"begA1"\t"endA1;
 
 
 	### 2) Block A2

 	## Block beginning
 	# A) 5-mer beginning higher than chromosome end coord (Upper bound)
 	if ((endA+1) > chrLengths[chrA])
 	{
 		begA2=chrLengths[chrA];
 	}
 	# B) 5-mer beginning within chromosomal coord.
 	else
 	{
 		begA2=(endA+1);
 	}
 	
 	## Block end

 	# A) 5-mer beginning higher than chromosome end coord (Upper bound)
 	if ((endA+5) > chrLengths[chrA])
 	{
 		endA2=chrLengths[chrA];
 	}
 	# C) 5-mer beginning within chromosomal coord.
 	else
 	{
 		endA2=(endA+5);
 	}
 	
 	## Block coordinates 	
 	blockA2=chrA"\t"strA"\t"begA2"\t"endA2;
 	
 	### 3) Block B1
 	
 	## Block beginning	
 	# A) 5-mer beginning lower than chromosome beginning coord (lower bound)
 	if ((begB-5) < 1)
 	{
 		begB1="1";
 	}
 	# B) 5-mer beginning within chromosomal coord.
 	else
 	{
 		begB1=(begB-5);
 	}
 	
	## Block end
 	# A) 5-mer beginning lower than chromosome beginning coord (lower bound)
 	if ((begB-1) < 1)
 	{
 		endB1="1";
 	}
 	# B) 5-mer beginning within chromosomal coord.
 	else
 	{
 		endB1=(begB-1);
 	}
 	
 	blockB1=chrB"\t"strB"\t"begB1"\t"endB1;
 	
	### 4) Block B2
	
 	## Block beginning	
 	# A) 5-mer beginning higher than chromosome end coord (Upper bound)
 	if ((endB+1) > chrLengths[chrB])
 	{
 		begB2=chrLengths[chrB];
 	}
 	# B) 5-mer beginning within chromosomal coord.
 	else
 	{
 		begB2=(endB+1);
 	}
 	
 	## Block end	
 	# B) 5-mer beginning higher than chromosome end coord (Upper bound)
 	if ((endB+5) > chrLengths[chrB])
 	{
 		endB2=chrLengths[chrB];
 	}
 	# C) 5-mer beginning within chromosomal coord.
 	else
 	{
 		endB2=(endB+5);
 	}
 	
 	blockB2=chrB"\t"strB"\t"begB2"\t"endB2;
 	
 	#### Print block coordinates
 	print	blockA1;
 	print	blockA2;
 	print	blockB1;
 	print	blockB2;
    
}

