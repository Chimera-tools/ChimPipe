#!/usr/bin/env awk

# *****************************************************************************
	
#	ChimFilter.awk
	
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
# Filters chimeric junction candidates to produce a final set of chimeric junctions. 

# Usage
#######

# awk -f chimjunc_filter.awk chimeric_junctions_candidates.txt

# Input:
# juncCoord	totalNbPE	nbSpanningPE	nbStag	percStag	nbMulti	percMulti	nbDiscordantPE	nbInconsistentPE	percInconsistentPE	overlapA	overlapB	distExonBoundaryA	distExonBoundaryB	maxBLastAlignLen	blastAlignSim	donorSS	acceptorSS	beg	end	sameChrStr	okGxOrder	dist	gnIdA	 gnIdB	gnNameA	gnNameB	gnTypeA	gnTypeB	juncSpanningReadIds	supportingPairsIds	inconsistentPairsIds
#chr21_42880008_-:chr21_39817544_-	200	126	32	25.3968	36	28.5714	74	14	15.9091	100	100	0	0	na	na	GT	AG	42880059	39817494	1	1	3062464	ENSG00000184012.7	ENSG00000157554.14	TMPRSS2	ERG	protein_coding	protein_coding

BEGIN{
	
	## Set default filtering parameters configuration in case they do not have been set by the user
	
	# 1. Minimum number of total paired-ends (staggered + discordant)
	
	if (minTotalNbPE == "")
	{
		minTotalNbPE=3;
	}
	
	# 2. Minimum number of spanning paired-ends 
	
	if (minNbSpanningPE == "")
	{
		minNbSpanningPE=1;
	}
	
	# 3. Minimum number of discordant paired-ends 
	
	if (minNbDiscordantPE == "")
	{
		minNbDiscordantPE=1;
	}
	
	# 4. Minimum percentage of staggered reads
	
	if (minPercStaggered == "")
	{
		minPercStaggered=0;
	}
	
	# 5. Maximum percentage of multimapped spanning reads
	
	if (maxPercMultimaps == "")
	{
		maxPercMultimaps=95;
	}
	
	# 6. Maximum percentage of inconsistent paired ends
	
	if (maxPercInconsistentPE == "")
	{
		maxPercInconsistentPE=50;
	}
	
	# 7. Maximum distance to exon boundary
	
	if (maxDistExonBoundary == "")
	{
		maxDistExonBoundary="-1";
	}
	
	## Print header:
	header="juncCoord\tfiltered\treason\ttotalNbPE\tnbSpanningPE\tnbStag\tpercStag\tnbMulti\tpercMulti\tnbDiscordantPE\tnbInconsistentPE\tpercInconsistentPE\toverlapA\toverlapB\tdistExonBoundaryA\tdistExonBoundaryB\tmaxBLastAlignLen\tblastAlignSim\tdonorSS\tacceptorSS\tbeg\tend\tsameChrStr\tokGxOrder\tdist\tgnIdA\tgnIdB\tgnNameA\tgnNameB\tgnTypeA\tgnTypeB\tjuncSpanningReadIds\tsupportingPairsIds\tinconsistentPairsIds";
	
	print header;
}

NR>1{
	# Get input information
	juncCoord=$1;	
	totalNbPE=$2;	
	nbSpanningPE=$3;	
	nbStag=$4;	
	percStaggered=$5;	
	nbMultimaps=$6;	
	percMultimaps=$7;	
	nbDiscordantPE=$8;	
	nbInconsistentPE=$9;	
	percInconsistentPE=$10;	
	overlapA=$11;	
	overlapB=$12;	
	distExonBoundaryA=$13;	
	distExonBoundaryB=$14;	
	maxBLastAlignLen=$15;	
	blastAlignSim=$16;
	donorSS=$17;	
	acceptorSS=$18;	
	beg=$19;	
	end=$20;	
	sameChrStr=$21;	
	okGxOrder=$22;	
	dist=$23;	
	gnIdA=$24;	
	gnIdB=$25;	
	gnNameA=$26;	
	gnNameB=$27;	
	gnTypeA=$28;	
	gnTypeB=$29;	
	juncSpanningReadIds=$30;	
	supportingPairsIds=$31;	
	inconsistentPairsIds=$32;
	
	filtered=0;
	reason="";
	
	## Filter 1: Minimum number of total paired-ends (staggered + discordant)
	if (totalNbPE <= minTotalNbPE)
	{
		filtered=1;
		reason="totalPE,"reason
	}
	
	## Filter 2: Minimum number of spanning paired-ends 
	if (nbSpanningPE <= minNbSpanningPE)
	{
		filtered=1;
		reason="spanningPE,"reason
	}
	
	## Filter 3: Minimum number of discordant paired-ends 
	if (nbDiscordantPE <= minNbDiscordantPE)
	{
		filtered=1;
		reason="discordantPE,"reason
	}
	
	## Filter 4: Minimum percentage of staggered reads
	if (percStaggered <= minPercStaggered) 
	{
		filtered=1;
		reason="percStag,"reason
	}	
	
	## Filter 5: Maximum percentage of multimapped spanning reads
	if (percMultimaps > maxPercMultimaps)
	{
		filtered=1;
		reason="percMulti,"reason
	}
	
	## Filter 6: Maximum percentage of inconsistent paired ends
	if ((percInconsistentPE > maxPercInconsistentPE) && (percInconsistentPE != "na"))
	{
		filtered=1;
		reason="percInconsistentPE,"reason
	}
	
	## Filter 7: Maximum distance to exon boundary
	if (((distExonBoundaryA > maxDistExonBoundary) || (distExonBoundaryB > maxDistExonBoundary)) && (maxDistExonBoundary != "-1"))
	{
		filtered=1;
		reason="distExonBoundary,"reason
	}
	
	if (reason == "")
	{
		reason="na";
	}
	
	## Print output
row=juncCoord"\t"filtered"\t"reason"\t"totalNbPE"\t"nbSpanningPE"\t"nbStag"\t"percStaggered"\t"nbMultimaps"\t"percMultimaps"\t"nbDiscordantPE"\t"nbInconsistentPE"\t"percInconsistentPE"\t"overlapA"\t"overlapB"\t"distExonBoundaryA"\t"distExonBoundaryB"\t"maxBLastAlignLen"\t"blastAlignSim"\t"donorSS"\t"acceptorSS"\t"beg"\t"end"\t"sameChrStr"\t"okGxOrder"\t"dist"\t"gnIdA"\t"gnIdB"\t"gnNameA"\t"gnNameB"\t"gnTypeA"\t"gnTypeB"\t"juncSpanningReadIds"\t"supportingPairsIds"\t"inconsistentPairsIds;		
	print row; 

}
