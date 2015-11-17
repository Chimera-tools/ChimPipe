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
	
	#### Set default filtering parameters configuration in case they do not have been set by the user
	
	### 1. Filters based on the number of supporting reads
	# A) Junction in exon boundaries (annotated splice-sites)
	# A.1. Minimum number of total supporting evidences (spanning reads + consistent PE)
	
	if (minNbTotal == "")
	{
		minNbTotal=4;
	}
	
	# A.2. Minimum number of spanning reads 
	if (minNbSpanning == "")
	{
		minNbSpanning=1;
	}
	
	# A.3. Minimum number of consistent paired-ends 
	if (minNbConsistentPE == "")
	{
		minNbConsistentPE=2;
	}
	
	# B) Junction not in exon boundaries (at least one novel splice-sites)
	# B.1. Minimum number of total supporting evidences (spanning reads + consistent PE)
	if (minNbTotalNovelSS == "")
	{
		minNbTotalNovelSS=8;
	}
	
	# B.2. Minimum number of spanning reads 
	if (minNbSpanningNovelSS == "")
	{
		minNbSpanningNovelSS=2;
	}
	
	# B.3. Minimum number of consistent paired-ends 
	if (minNbConsistentPENovelSS == "")
	{
		minNbConsistentPENovelSS=4;
	}
	
	### 2. Filters based on percentages
	# 2.1. Minimum percentage of staggered reads
	if (minPercStaggered == "")
	{
		minPercStaggered=0;
	}
	
	# 2.2. Maximum percentage of multimapped spanning reads
	if (maxPercMultimaps == "")
	{
		maxPercMultimaps=100;
	}
	
	# 2.3. Maximum percentage of inconsistent paired ends
	if (maxPercInconsistentPE == "")
	{
		maxPercInconsistentPE=100;
	}
	
	### 3. Similarity based filter
	# A) Filtering configuration not provided
	if (similarityConf == "")
	{
		maxAlignLen=30;
		maxSimPerc=90;
	}
	# B) Filtering configuration provided
	else
	{
		split(similarityConf, similarityConfList, "+");
		maxAlignLen=similarityConfList[1];
		maxSimPerc=similarityConfList[2];
	} 
	
		
	## Print header:
	header="juncCoord\tfiltered\treason\tnbTotal(spanning+consistent)\tnbSpanningReads\tnbStaggered\tpercStaggered\tnbMulti\tpercMulti\tnbConsistentPE\tnbInconsistentPE\tpercInconsistentPE\toverlapA\toverlapB\tdistExonBoundaryA\tdistExonBoundaryB\tblastAlignLen\tblastAlignSim\tdonorSS\tacceptorSS\tbeg\tend\tsameChrStr\tokGxOrder\tdist\tgnIdsA\tgnIdsB\tgnNamesA\tgnNamesB\tgnTypesA\tgnTypesB\tjuncSpanningReadsIds\tconsistentPEIds\tinconsistentPEIds";
	
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
	
	filtered="0";
	reason="";
	
	
	# Apply filters
	################
	### 1. Filters based on the number of supporting reads
	## A) Junction in exon boundaries (annotated splice-sites)
	if ((distExonBoundaryA == "0") && (distExonBoundaryB == "0"))  
	{		
		# A.1. Minimum number of total supporting evidences (spanning reads + consistent PE)
		if (nbTotal < minNbTotal)
		{
			reason="totalSupport,"reason
		}
	
		# A.2. Minimum number of spanning reads
		if (nbSpanning < minNbSpanning)
		{
			reason="spanningReads,"reason
		}
	
		# A.3. Minimum number of consistent paired-ends 
		if (nbConsistentPE < minNbConsistentPE)
		{
			reason="consistentPE,"reason
		}
	}	
	## B) Junction not in exon boundaries (at least one novel splice-sites)
	else
	{	
		# B.1. Minimum number of total supporting evidences (spanning reads + consistent PE)
		if (nbTotal < minNbTotalNovelSS)
		{
			reason="totalSupport,"reason
		}
	
		# B.2. Minimum number of spanning reads
		if (nbSpanning < minNbSpanningNovelSS)
		{
			reason="spanningReads,"reason
		}
		
		# B.3. Minimum number of consistent paired-ends 
		if (nbConsistentPE < minNbConsistentPENovelSS)
		{
			reason="consistentPE,"reason
		}
	}
	
	### 2. Filters based on percentages
	## 2.1 Minimum percentage of staggered reads
	if (percStaggered < minPercStaggered) 
	{
		reason="percStag,"reason
	}	
	
	## 2.2 Maximum percentage of multimapped spanning reads
	if (percMultimaps > maxPercMultimaps)
	{
		reason="percMulti,"reason
	}
	
	## 2.3 Maximum percentage of inconsistent paired ends
	if ((percInconsistentPE > maxPercInconsistentPE) && (percInconsistentPE != "na"))
	{
		reason="percInconsistentPE,"reason
	}
	
	### 3. Similarity based filter	
	if ((blastAlignLen > maxAlignLen) && (blastAlignSim > maxSimPerc) && (blastAlignLen != "na"))
	{
		reason="similarity,"reason
	}

	# Set flags according filtering outcome
	########################################
	
	# A) Junction pass all the filters
	if (reason == "")
	{
		reason="na";
	}
	# B) Junction do not pass any filter
	else
	{
		filtered="1";
	}
	
	# Print output
	###############
row=juncCoord"\t"filtered"\t"reason"\t"nbTotal"\t"nbSpanning"\t"nbStaggered"\t"percStaggered"\t"nbMultimaps"\t"percMultimaps"\t"nbConsistentPE"\t"nbInconsistentPE"\t"percInconsistentPE"\t"overlapA"\t"overlapB"\t"distExonBoundaryA"\t"distExonBoundaryB"\t"blastAlignLen"\t"blastAlignSim"\t"donorSS"\t"acceptorSS"\t"beg"\t"end"\t"sameChrStr"\t"okGxOrder"\t"dist"\t"gnIdsA"\t"gnIdsB"\t"gnNamesA"\t"gnNamesB"\t"gnTypesA"\t"gnTypesB"\t"juncSpanningReadsIds"\t"consistentPEIds"\t"inconsistentPEIds;		
	print row; 

}
