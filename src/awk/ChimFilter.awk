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
# juncCoord	type	nbTotal(spanning+consistent)	nbSpanningReads	nbStaggered	percStaggered	nbMulti	percMulti	nbConsistentPE	nbInconsistentPE	percInconsistentPE	overlapA	overlapB	distExonBoundaryA	distExonBoundaryB	blastAlignLen	blastAlignSim	donorSS	acceptorSS	beg	end	sameChrStr	okGxOrder	dist	gnIdsA	gnIdsB	gnNamesA	gnNamesB	gnTypesA	gnTypesB	juncSpanningReadsIds	consistentPEIds	inconsistentPEIds
# chr19_57176129_+:chr19_58879504_+	intrachromosomal	2	2	2	100	2	100	0	0	na	100	100484	514	133	81.95	GC	AG	57176088	58879518	1	1	1703375	ENSG00000127903.12	ENSG00000152475.6	ZNF835	ZNF837	protein_coding	protein_coding	SRR201779.989833_PATHBIO-SOLEXA2_30TUEAAXX:3:15:325:1075_length=53#0/1,SRR201779.3705596_PATHBIO-SOLEXA2_30TUEAAXX:3:57:1642:1097_length=53#0/2,	na	na


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
	
	### 4. Biotype blacklist filter
	if (biotype == "")
	{
		biotype="pseudogene";
	}
		
	## Print header:
	header="juncCoord\ttype\tfiltered\treason\tnbTotal(spanning+consistent)\tnbSpanningReads\tnbStaggered\tpercStaggered\tnbMulti\tpercMulti\tnbConsistentPE\tnbInconsistentPE\tpercInconsistentPE\toverlapA\toverlapB\tdistExonBoundaryA\tdistExonBoundaryB\tblastAlignLen\tblastAlignSim\tdonorSS\tacceptorSS\tbeg\tend\tsameChrStr\tokGxOrder\tdist\tgnIdsA\tgnIdsB\tgnNamesA\tgnNamesB\tgnTypesA\tgnTypesB\tjuncSpanningReadsIds\tconsistentPEIds\tinconsistentPEIds";
	
	print header;
}


NR>1{	
	
	# Get input information
	########################
	juncCoord=$1;	
	type=$2;
	nbTotal=$3;	
	nbSpanning=$4;	
	nbStaggered=$5;	
	percStaggered=$6;	
	nbMultimaps=$7;	
	percMultimaps=$8;	
	nbConsistentPE=$9;	
	nbInconsistentPE=$10;	
	percInconsistentPE=$11;	
	overlapA=$12;	
	overlapB=$13;	
	distExonBoundaryA=$14;	
	distExonBoundaryB=$15;	
	blastAlignLen=$16;	
	blastAlignSim=$17;
	donorSS=$18;	
	acceptorSS=$19;	
	beg=$20;	
	end=$21;	
	sameChrStr=$22;	
	okGxOrder=$23;	
	dist=$24;	
	gnIdsA=$25;	
	gnIdsB=$26;	
	gnNamesA=$27;	
	gnNamesB=$28;	
	gnTypesA=$29;	
	gnTypesB=$30;	
	juncSpanningReadsIds=$31;	
	consistentPEIds=$32;	
	inconsistentPEIds=$33;
	
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

	### 4. Biotype blacklist filter
	
	nb=split(biotype, biotypeBlackList, ",");
	
	nbA=split(gnTypesA, biotypesListA, ",");
	nbB=split(gnTypesB, biotypesListB, ",");
	
	
	## A) Iterate over the biotype of the genes with exons overlapping the junction donor side
	for (i=1;i<=nbA;i++)
	{
		if (biotypesListA[i] != "")
		{
			# Iterate over the biotypes in the black list
			for (j=1;j<=nb;j++)
			{
				# Discard the junction if it involves genes with a biotypes in the black list
				if (biotypesListA[i] == biotypeBlackList[j])
				{
					reason="biotype,"reason
				}
			}	
		}
	}
	
	## B) Iterate over the biotype of the genes with exons overlapping the junction acceptor side
	for (i=1;i<=nbB;i++)
	{
		if (biotypesListB[i] != "")
		{
			# Iterate over the biotypes in the black list
			for (j=1;j<=nb;j++)
			{
				# Discard the junction if it involves genes with a biotypes in the black list
				if (biotypesListB[i] == biotypeBlackList[j])
				{
					reason="biotype,"reason
				}
			}	
		}
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
row=juncCoord"\t"type"\t"filtered"\t"reason"\t"nbTotal"\t"nbSpanning"\t"nbStaggered"\t"percStaggered"\t"nbMultimaps"\t"percMultimaps"\t"nbConsistentPE"\t"nbInconsistentPE"\t"percInconsistentPE"\t"overlapA"\t"overlapB"\t"distExonBoundaryA"\t"distExonBoundaryB"\t"blastAlignLen"\t"blastAlignSim"\t"donorSS"\t"acceptorSS"\t"beg"\t"end"\t"sameChrStr"\t"okGxOrder"\t"dist"\t"gnIdsA"\t"gnIdsB"\t"gnNamesA"\t"gnNamesB"\t"gnTypesA"\t"gnTypesB"\t"juncSpanningReadsIds"\t"consistentPEIds"\t"inconsistentPEIds;		
	print row; 

}
