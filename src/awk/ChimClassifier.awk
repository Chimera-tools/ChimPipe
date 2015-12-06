#!/usr/bin/env awk

# *****************************************************************************
	
#	ChimClassifier.awk
	
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
################ 
# Classify chimeric junctions according to the chromomal position of the donor and acceptor sites in one of these categories:

# 1) readthrough -	donor and acceptor sites within the same chromosome, strand and within less than X base pairs (150.000 default)

# 2) intrachromosomal - donor and acceptor sites within the same chromosome, strand and in separated by more than X base pairs (150.000 default)
			
# 3) inverted - donor site downstream than acceptor site (opposite as expected) and both in the same chromosome
 
# 4) interstrand - donor and acceptor sites within the same chromosome but in different strand

# 5) interchromosomal -	donor and acceptor sites in different chromosomes
		

# Usage
#######

# awk -f ChimClassifier.awk chimeric_junctions_candidates.txt

# Input:
# juncCoord	totalNbPE	nbSpanningPE	nbStag	percStag	nbMulti	percMulti	nbDiscordantPE	nbInconsistentPE	percInconsistentPE	overlapA	overlapB	distExonBoundaryA	distExonBoundaryB	maxBLastAlignLen	blastAlignSim	donorSS	acceptorSS	beg	end	sameChrStr	okGxOrder	dist	gnIdA	 gnIdB	gnNameA	gnNameB	gnTypeA	gnTypeB	juncSpanningReadIds	supportingPairsIds	inconsistentPairsIds
#chr21_42880008_-:chr21_39817544_-	200	126	32	25.3968	36	28.5714	74	14	15.9091	100	100	0	0	na	na	GT	AG	42880059	39817494	1	1	3062464	ENSG00000184012.7	ENSG00000157554.14	TMPRSS2	ERG	protein_coding	protein_coding

BEGIN{
	
	#### Set default classifier configuration
	
	if (readthroughMaxDist == "")
	{
		readthroughMaxDist=150000;
	}
		
	## Print header:
	header="juncCoord\ttype\tnbTotal(spanning+consistent)\tnbSpanningReads\tnbStaggered\tpercStaggered\tnbMulti\tpercMulti\tnbConsistentPE\tnbInconsistentPE\tpercInconsistentPE\toverlapA\toverlapB\tdistExonBoundaryA\tdistExonBoundaryB\tblastAlignLen\tblastAlignSim\tdonorSS\tacceptorSS\tbeg\tend\tsameChrStr\tokGxOrder\tdist\tgnIdsA\tgnIdsB\tgnNamesA\tgnNamesB\tgnTypesA\tgnTypesB\tjuncSpanningReadsIds\tconsistentPEIds\tinconsistentPEIds";
	
	print header;
}

NR>1{
	
	# Get input information
	########################
	juncCoord=$1;	
	nbTotal=$2;	
	nbSpanning=$3;	
	nbStaggered=$4;	
	percStaggered=$5;	
	nbMultimaps=$6;	
	percMultimaps=$7;	
	nbConsistentPE=$8;	
	nbInconsistentPE=$9;	
	percInconsistentPE=$10;	
	overlapA=$11;	
	overlapB=$12;	
	distExonBoundaryA=$13;	
	distExonBoundaryB=$14;	
	blastAlignLen=$15;	
	blastAlignSim=$16;
	donorSS=$17;	
	acceptorSS=$18;	
	beg=$19;	
	end=$20;	
	sameChrStr=$21;	
	okGxOrder=$22;	
	dist=$23;	
	gnIdsA=$24;	
	gnIdsB=$25;	
	gnNamesA=$26;	
	gnNamesB=$27;	
	gnTypesA=$28;	
	gnTypesB=$29;	
	juncSpanningReadsIds=$30;	
	consistentPEIds=$31;	
	inconsistentPEIds=$32;
		
			
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
row=juncCoord"\t"type"\t"nbTotal"\t"nbSpanning"\t"nbStaggered"\t"percStaggered"\t"nbMultimaps"\t"percMultimaps"\t"nbConsistentPE"\t"nbInconsistentPE"\t"percInconsistentPE"\t"overlapA"\t"overlapB"\t"distExonBoundaryA"\t"distExonBoundaryB"\t"blastAlignLen"\t"blastAlignSim"\t"donorSS"\t"acceptorSS"\t"beg"\t"end"\t"sameChrStr"\t"okGxOrder"\t"dist"\t"gnIdsA"\t"gnIdsB"\t"gnNamesA"\t"gnNamesB"\t"gnTypesA"\t"gnTypesB"\t"juncSpanningReadsIds"\t"consistentPEIds"\t"inconsistentPEIds;		
	print row; 

}
