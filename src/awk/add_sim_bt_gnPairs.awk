#!/usr/bin/env awk

# *****************************************************************************
	
#	add_sim_bt_gnPairs.awk
	
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
# Need to be done..

## Input:

# 1) chimeric_spliceJunctions_discordantPE.txt
# juncCoord	totalNbPE	nbSpanningPE	nbStag	percStag	nbMulti	percMulti	nbDiscordantPE	nbInconsistentPE	percInconsistentPE	overlapA	overlapB	distExonBoundaryA	distExonBoundaryB	donorSS	acceptorSS	beg	end	sameChrStr	okGxOrder	dist	gnIdA	 gnIdB	gnNameA	gnNameB	gnTypeA	gnTypeB	juncSpanningReadIds	supportingPairsIds	inconsistentPairsIds
# chr21_42880008_-:chr21_39817544_-	200	126	32	25.3968	36	28.5714	74	14	15.9091	100	100	0	0	GT	AG	42880059	39817494	1	1	3062464	ENSG00000184012.7	ENSG00000157554.14	TMPRSS2	ERG	protein_coding	protein_coding ...

# 2) fileRef=annotId_gene1_gene2_alphaorder_pcentsim_lgalign_trpair.txt
# ENSG00000000003.10 ENSG00000042317.12 100.00 31 ENST00000373020.4 ENST00000553885.1
# ENSG00000000003.10 ENSG00000053524.7 91.67 84 ENST00000373020.4 ENST00000461074.1
# ENSG00000000003.10 ENSG00000059804.11 87.50 80 ENST00000373020.4 ENST00000075120.7

BEGIN{
	while (getline < fileRef >0)
	{
		# Declare variables with relevant input information
		gnIdA=$1;
		gnIdB=$2;
		percSim=$3;
		alignLength=$4;
		
		alignLengths[gnIdA"-"gnIdB]=alignLength;
		percSimilarities[gnIdA"-"gnIdB]=percSim;
	}
}

# Print header
NR==1{
	header="juncCoord\ttotalNbPE\tnbSpanningPE\tnbStag\tpercStag\tnbMulti\tpercMulti\tnbDiscordantPE\tnbInconsistentPE\tpercInconsistentPE\toverlapA\toverlapB\tdistExonBoundaryA\tdistExonBoundaryB\tmaxBLastAlignLen\tblastAlignSim\tdonorSS\tacceptorSS\tbeg\tend\tsameChrStr\tokGxOrder\tdist\tgnIdA\t gnIdB\tgnNameA\tgnNameB\tgnTypeA\tgnTypeB\tjuncSpanningReadIds\tsupportingPairsIds\tinconsistentPairsIds";
    print header;
}

NR>1{
	# Declare variables with chimeric junctions plus associated information
	juncCoord=$1;
	totalNbPE=$2;
	nbSpanningPE=$3;
	nbStag=$4;
	percStag=$5;
	nbMulti=$6;
	percMulti=$7;
	nbDiscordantPE=$8;
	nbInconsistentPE=$9;
	percInconsistentPE=$10;
	overlapA=$11;
	overlapB=$12;
	distExonBoundaryA=$13;
	distExonBoundaryB=$14;
	donorSS=$15;
	acceptorSS=$16;
	beg=$17;
	end=$18;
	sameChrStr=$19;
	okGxOrder=$20;
	dist=$21;
	gnIdA=$22;
	gnIdB=$23;
	gnNameA=$24;
	gnNameB=$25;
	gnTypeA=$26;
	gnTypeB=$27;
	juncSpanningReadIds=$28;
	supportingPairsIds=$29;
	inconsistentPairsIds=$30;
	
	maxBLastAlignLen=0;
	blastAlignSim=0; 	
	
	##
	split(gnIdA, gnlistA, ",");
	split(gnIdB,gnlistB,",");
	
	for (gnA in gnlistA)
	{
		for (gnB in gnlistB)
		{
			if ((gnlistA[gnA]!="") || (gnlistB[gnB]!=""))
			{
				gnPairIds=((gnlistA[gnA]<=gnlistB[gnB]) ? (gnlistA[gnA]"-"gnlistB[gnB]) : (gnlistB[gnB]"-"gnlistA[gnA]));
				
				if ((alignLengths[gnPairIds] > maxBLastAlignLen) || ((alignLengths[gnPairIds] == maxBLastAlignLen) && (percSimilarities[gnPairIds] > blastAlignSim)))
				{
					maxBLastAlignLen=alignLengths[gnPairIds];
					blastAlignSim=percSimilarities[gnPairIds];
				}
			}
		}
	}
	
	if ((maxBLastAlignLen == 0) || (blastAlignSim == 0))	 
	{
		maxBLastAlignLen="na";
		blastAlignSim="na";
	}
	
 
 	## Print output
 	
		  row=juncCoord"\t"totalNbPE"\t"nbSpanningPE"\t"nbStag"\t"percStag"\t"nbMulti"\t"percMulti"\t"nbDiscordantPE"\t"nbInconsistentPE"\t"percInconsistentPE"\t"overlapA"\t"overlapB"\t"distExonBoundaryA"\t"distExonBoundaryB"\t"maxBLastAlignLen"\t"blastAlignSim"\t"donorSS"\t"acceptorSS"\t"beg"\t"end"\t"sameChrStr"\t"okGxOrder"\t"dist"\t"gnIdA"\t"gnIdB"\t"gnNameA"\t"gnNameB"\t"gnTypeA"\t"gnTypeB"\t"juncSpanningReadIds"\t"supportingPairsIds"\t"inconsistentPairsIds;
	        
	print row;
}



