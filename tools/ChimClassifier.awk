#!/usr/bin/env awk

# *****************************************************************************
	
#	ChimClassifier.awk
	
#	This file is part of the ChimPipe pipeline 

#	Copyright (c) 2014-2016            Bernardo Rodríguez-Martín 
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
################ 
# Classify chimeric junctions according to the chromomal position of the donor and acceptor sites in one of these categories:

# 1) readthrough - donor and acceptor sites within the same chromosome, strand and within less than X base pairs (150.000 default)
# 2) intrachromosomal - donor and acceptor sites within the same chromosome, strand and in separated by more than X base pairs (150.000 default)
# 3) inverted - donor site > acceptor site coordinates (opposite as expected) and both in the same chromosome
# 4) interstrand - donor and acceptor sites within the same chromosome but in different strand
# 5) interchromosomal -	donor and acceptor sites in different chromosomes

# Usage
#######

# awk -f ChimClassifier.awk chimJunctions.txt

### Input:
## txt with chimeric junctions in ChimPipe format (see 'juncCoord' at https://chimpipe.readthedocs.io/en/latest/manual.html#output)
# chrY_14802370_+:chrY_14813939_+       
# chr10_51387763_+:chr10_51732772_+     
# chr5_170720985_+:chr5_169267761_+     
# chr2_99193606_+:chr2_234746298_-      
# chr16_85023910_-:chr12_123444869_-    

### Output:
## txt with chimeric junctions classified 
# chrY_14802370_+:chrY_14813939_+       readthrough
# chr10_51387763_+:chr10_51732772_+	intrachromosomal
# chr5_170720985_+:chr5_169267761_+     inverted
# chr2_99193606_+:chr2_234746298_-	interstrand
# chr16_85023910_-:chr12_123444869_-    interchromosomal

BEGIN{
	
	#### Set default classifier configuration
	
	if (readthroughMaxDist == "")
	{
		readthroughMaxDist=150000;
	}
		
	## Print header:
	header="juncCoord\ttype";	
	print header;
}

{
	
	# Get input information
	########################
	juncCoord=$1;	
			
	# Classify chimeric junctions
	###############################
	
	split(juncCoord,a,":");
	split(a[1],b1,"_");
	split(a[2],b2,"_");
	
	chrA=b1[1];
	coordA=b1[2];
	strandA=b1[3];
	chrB=b2[1];
	coordB=b2[2];
	strandB=b2[3]; 

	### A) Same chromosomes and strands
	if ((chrA == chrB) && (strandA == strandB))
	{
		### A.A) Expected genomic order
		if (((strandA == "+") && (coordA < coordB)) || ((strandA == "-") && (coordA > coordB)))
		{
			dist = coordB - coordA
			
			## A.A.A) Readthrough
			if (dist <= readthroughMaxDist)
			{
				type="readthrough"	
			}
			## A.A.B) Not readthrough
			else
			{
				type="intrachromosomal"
			}
		}
		### A.B) Unexpected genomic order
		else
		{	
		type="inverted"
		}
	}
	### B) Same chromosomes and different strands
	else if ((chrA == chrB) && (strandA != strandB))
	{
		type="interstrand"
	}
	### C) Different chromosomes
	else
	{
		type="interchromosomal"
	}
	
	# Print output
	###############
	row=juncCoord"\t"type;		
	print row; 

}
