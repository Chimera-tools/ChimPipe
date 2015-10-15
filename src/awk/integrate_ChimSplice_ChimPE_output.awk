#!/usr/bin/env awk

# *****************************************************************************

#       integrate_ChimSplice_ChimPE_output.awk

#       This file is part of the ChimPipe pipeline 

#       Copyright (c) 2014 Bernardo Rodríguez-Martín 
#                                          Emilio Palumbo 
#                                          Sarah djebali 

#       Computational Biology of RNA Processing group
#       Department of Bioinformatics and Genomics
#       Centre for Genomic Regulation (CRG)
                                           
#       Github repository - https://github.com/Chimera-tools/ChimPipe

#       Documentation - https://chimpipe.readthedocs.org/

#       Contact - chimpipe.pipeline@gmail.com

#       Licenced under the GNU General Public License 3.0 license.
#******************************************************************************

# Description
###############

# Inputs:
##########
## 1) ChimSplice output file: chimeric_spliceJunctions.txt
# chr21_42880008_-:chr21_39817544_-       CHIMERIC        126     32      25.3968 36      28.5714 100     100     0       0       GT      AG      42880059        39817494        1       1   3062464  ENSG00000184012.7       ENSG00000157554.14      TMPRSS2 ERG     protein_coding  protein_coding

## 2) ChimPE output file: discordant_readPairs.txt
# SRR201779.6360963_PATHBIO-SOLEXA2_30TUEAAXX:3:98:699:447_length=53#0    DISCORDANT      chr21_39817369_39817421_+       chr21_42880008_42880060_-       ENSG00000157554.14      ENSG00000184012.7
# SRR201779.6043897_PATHBIO-SOLEXA2_30TUEAAXX:3:93:627:338_length=53#0    DISCORDANT      chr21_39817370_39817422_+       chr21_42880014_42880066_-       ENSG00000157554.14      ENSG00000184012.7
# SRR201779.1492533_PATHBIO-SOLEXA2_30TUEAAXX:3:22:1310:1233_length=53#0  DISCORDANT      chr21_39817374_39817426_+       chr21_42880014_42880066_-       ENSG00000157554.14      ENSG00000184012.7

## output:


# Usage example:
################
# awk -v ChimPE=chimJunctions_chimSplice.txt -f integrate_ChimSplice_ChimPE_output.awk  discordantPE.txt


function discordantPEsupport(pairId, donorGene, acceptorGene, mapCoordMateSupDonor, mapCoordMateSupAcceptor, chimJunctions)
{
	# 1) Get mapping coordinates of discordant paired end mates
	# 1.1) Mate supporting donor
	nbBlocksDonor=split(mapCoordMateSupDonor, mapblocks, ":")
        			
    # 1.1.A) Mate mapped in one block
    if (nbBlocksDonor == "1")
    {
     	split(mapCoordMateSupDonor, coords, "_");
	    begMateSupDonor=coords[2];
		endMateSupDonor=coords[3];	
	}
	# 1.1.B) Split-mapped mate
    else
    {
    	split(mapblocks[1], coordsA, "_");
    	split(mapblocks[2], coordsB, "_");
    				
    	strand=coordsA[3];
    				
    	# splice junction between exons in plus strand
    	if (strand=="+")
    	{
    		begMateSupDonor=coordsA[2];
    		endMateSupDonor=coordsB[2];
    	}
    	# splice junction between exons in minus strand
    	else
    	{
    		begMateSupDonor=coordsB[2];
    		endMateSupDonor=coordsA[2];
    	}
    }    	
		
	# 1.2) Mate supporting acceptor
	nbBlocksDonor=split(mapCoordMateSupAcceptor, mapblocks, ":")
        			
    # 1.2.A) Mate mapped in one block
    if (nbBlocksDonor == "1")
    {
    	split(mapCoordMateSupAcceptor, coords, "_");
		begMateSupAcceptor=coords[2];
		endMateSupAcceptor=coords[3];	
	}
	# 1.2.B) Split-mapped mate
    else
    {
    	split(mapblocks[1], coordsA, "_");
    	split(mapblocks[2], coordsB, "_");
    				
    	strand=coordsA[3];
    				
    	# splice junction between exons in plus strand
    	if (strand=="+")
    	{
    		begMateSupAcceptor=coordsA[2];
    		endMateSupAcceptor=coordsB[2];
    	}
    	# splice junction between exons in minus strand
    	else
    	{
    		begMateSupAcceptor=coordsB[2];
    		endMateSupAcceptor=coordsA[2];
    	}
    }    				        		
        				
    #print $0;
    #print "juncList:" chimJunctions, donorGene, acceptorGene;
    #print chimJunctions, begMateSupDonor, endMateSupDonor, begMateSupAcceptor, endMateSupAcceptor;
        		
        		
    #print chimJunctionsList[gnListA[i]"-"gnListB[j]];
   
   
   	# 2) Iterate over the list of chimeric splice junctions connecting GnA-GnB or GnB-GnA
   	# For each splice junction check if the discordant paired-end map trully flank the junction or not       		
   	nbChimJunc=split(chimJunctions, chimJuncList, ",")
        		
    for (k=1;k<=nbChimJunc;k++)
    {
     	if (chimJuncList[k] != "")
     	{
        	chimJunc=chimJuncList[k];
        				
        	split(chimJunc, chimJunctionSplit, ":");
        				
			split(chimJunctionSplit[1], donor, "_");
			split(chimJunctionSplit[2], acceptor,"_");
					
			coordDonor=donor[2];
			strandDonor=donor[3];
						
			coordAcceptor=acceptor[2];
			strandAcceptor=acceptor[3];
					
			# print chrDonor, coordDonor, strandDonor, chrAcceptor, coordAcceptor, strandAcceptor;
						
			chimJuncFlankingPE(pairId, coordDonor, strandDonor, coordAcceptor, strandAcceptor, begMateSupDonor, endMateSupDonor, begMateSupAcceptor, endMateSupAcceptor);					
		}
	}        		
	#print "########################";
}



function chimJuncFlankingPE(pairId, coordDonor, strandDonor, coordAcceptor, strandAcceptor, begMateSupDonor, endMateSupDonor, begMateSupAcceptor, endMateSupAcceptor)
{
	if (strandDonor=="+")
	{
		if (strandAcceptor=="+")  # Chimeric junction in +/+ 
		{
			## A) Discordant paired-end flanking chimeric splice junction
			if ((endMateSupDonor <= coordDonor) && (begMateSupAcceptor >= coordAcceptor))
			{
				nbSupportingPairs[donorGene"-"acceptorGene"::"chimJunc]++;
				supportingPairsIds[donorGene"-"acceptorGene"::"chimJunc]=pairId","supportingPairsIds[donorGene"-"acceptorGene"::"chimJunc];	
        				
				# print "junc_con", coordDonor, strandDonor, coordAcceptor, strandAcceptor;
        		# print "mates_con", begMateSupDonor, endMateSupDonor, begMateSupAcceptor, endMateSupAcceptor;
        		# print "nb supporting", nbSupportingPairs[donorGene"-"acceptorGene"::"chimJunc];
        	}
        	## B) Discordant paired-end do not flanking chimeric splice junction
        	else
        	{
        		nbInconsistentPairs[donorGene"-"acceptorGene"::"chimJunc]++;
        		inconsistentPairsIds[donorGene"-"acceptorGene"::"chimJunc]=pairId","inconsistentPairsIds[donorGene"-"acceptorGene"::"chimJunc];	
        				
        		#print "junc_incon", coordDonor, strandDonor, coordAcceptor, strandAcceptor;
        		#print "mates_incon", begMateSupDonor, endMateSupDonor, begMateSupAcceptor, endMateSupAcceptor;
        		#print "nb inconsistent", nbInconsistentPairs[donorGene"-"acceptorGene"::"chimJunc];
        	}
        }
        else  # Chimeric junction in +/-
		{
			## A) Discordant paired-end flanking chimeric splice junction
			if ((endMateSupDonor <= coordDonor) && (endMateSupAcceptor <= coordAcceptor))
			{
				nbSupportingPairs[donorGene"-"acceptorGene"::"chimJunc]++;
				supportingPairsIds[donorGene"-"acceptorGene"::"chimJunc]=pairId","supportingPairsIds[donorGene"-"acceptorGene"::"chimJunc];	
									
				#print "junc_con", coordDonor, strandDonor, coordAcceptor, strandAcceptor;
        		#print "mates_con", begMateSupDonor, endMateSupDonor, begMateSupAcceptor, endMateSupAcceptor;
        		#print "nb supporting", nbSupportingPairs[donorGene"-"acceptorGene"::"chimJunc];
        	}
        	## B) Discordant paired-end do not flanking chimeric splice junction
        	else
        	{
        		nbInconsistentPairs[donorGene"-"acceptorGene"::"chimJunc]++;
        		inconsistentPairsIds[donorGene"-"acceptorGene"::"chimJunc]=pairId","inconsistentPairsIds[donorGene"-"acceptorGene"::"chimJunc];	
        						
        		#print "junc_incon", coordDonor, strandDonor, coordAcceptor, strandAcceptor;
        		#print "mates_incon", begMateSupDonor, endMateSupDonor, begMateSupAcceptor, endMateSupAcceptor;
        		#print "nb inconsistent", nbInconsistentPairs[donorGene"-"acceptorGene"::"chimJunc];
        	}
		}
	}
	else
	{
		if (strandAcceptor=="+")  # Chimeric junction in -/+
		{
			## A) Discordant paired-end flanking chimeric splice junction
			if ((begMateSupDonor >= coordDonor) && (begMateSupAcceptor >= coordAcceptor))
			{
				nbSupportingPairs[donorGene"-"acceptorGene"::"chimJunc]++;		
				supportingPairsIds[donorGene"-"acceptorGene"::"chimJunc]=pairId","supportingPairsIds[donorGene"-"acceptorGene"::"chimJunc];	
									
				#print "junc_con", coordDonor, strandDonor, coordAcceptor, strandAcceptor;
    			#print "mates_con", begMateSupDonor, endMateSupDonor, begMateSupAcceptor, endMateSupAcceptor;
    			#print "nb supporting", nbSupportingPairs[donorGene"-"acceptorGene"::"chimJunc];
    		}
    		## B) Discordant paired-end do not flanking chimeric splice junction
        	else
        	{
        		nbInconsistentPairs[donorGene"-"acceptorGene"::"chimJunc]++;
        		inconsistentPairsIds[donorGene"-"acceptorGene"::"chimJunc]=pairId","inconsistentPairsIds[donorGene"-"acceptorGene"::"chimJunc];	
        							
        		#print "junc_incon", coordDonor, strandDonor, coordAcceptor, strandAcceptor;
        		#print "mates_incon", begMateSupDonor, endMateSupDonor, begMateSupAcceptor, endMateSupAcceptor;
        		#print "nb inconsistent", nbInconsistentPairs[donorGene"-"acceptorGene"::"chimJunc];
        	}
		}
		else  # Chimeric junction in -/-
		{	
			## A) Discordant paired-end flanking chimeric splice junction
			if ((begMateSupDonor >= coordDonor) && (endMateSupAcceptor <= coordAcceptor))
			{
				nbSupportingPairs[donorGene"-"acceptorGene"::"chimJunc]++;		
				supportingPairsIds[donorGene"-"acceptorGene"::"chimJunc]=pairId","supportingPairsIds[donorGene"-"acceptorGene"::"chimJunc];	
									
				#print "junc_con", coordDonor, strandDonor, coordAcceptor, strandAcceptor;
        		#print "mates_con", begMateSupDonor, endMateSupDonor, begMateSupAcceptor, endMateSupAcceptor;
        		#print "nb supporting", nbSupportingPairs[donorGene"-"acceptorGene"::"chimJunc];
        	}
        	## B) Discordant paired-end do not flanking chimeric splice junction
        	else
        	{
        		nbInconsistentPairs[donorGene"-"acceptorGene"::"chimJunc]++;
        		inconsistentPairsIds[donorGene"-"acceptorGene"::"chimJunc]=pairId","inconsistentPairsIds[donorGene"-"acceptorGene"::"chimJunc];	
        							
        		#print "junc_incon", coordDonor, strandDonor, coordAcceptor, strandAcceptor;
        		#print "mates_incon", begMateSupDonor, endMateSupDonor, begMateSupAcceptor, endMateSupAcceptor;
        		#print "nb inconsistent", nbInconsistentPairs[donorGene"-"acceptorGene"::"chimJunc];
        	}
		}
	}
}


BEGIN{
	header="juncCoord\ttotalNbPE\tnbSpanningPE\tnbStag\tpercStag\tnbMulti\tpercMulti\tnbDiscordantPE\tnbInconsistentPE\tpercInconsistentPE\toverlapA\toverlapB\tdistExonBoundaryA\tdistExonBoundaryB\tdonorSS\tacceptorSS\tbeg\tend\tsameChrStr\tokGxOrder\tdist\tgnIdA\t gnIdB\tgnNameA\tgnNameB\tgnTypeA\tgnTypeB\tjuncSpanningReadIds\tsupportingPairsIds\tinconsistentPairsIds";
	print header;
	
	
	while (getline < ChimSplice > 0)
	{	
		## Get all relevant information about a chimeric junction and save into variables or dictionaries
		
		#print $0;
		chimJunc=$1;
		chimJunctionsDic[chimJunc]="1";
		nbSpanningPE[chimJunc]=$3;
		nbStag[chimJunc]=$4;
		percStag[chimJunc]=$5;
		nbMulti[chimJunc]=$6;
		percMulti[chimJunc]=$7;
		overlapA[chimJunc]=$8;	
		overlapB[chimJunc]=$9;	
		distExonBoundaryA[chimJunc]=$10;	
		distExonBoundaryB[chimJunc]=$11;	
		donorSS[chimJunc]=$12;	
		acceptorSS[chimJunc]=$13;
		beg[chimJunc]=$14;	
		end[chimJunc]=$15;		
		sameChrStr[chimJunc]=$16;	
		okGxOrder[chimJunc]=$17;
		dist[chimJunc]=$18;	
		gnIdsA[chimJunc]=$19;	 
		gnIdsB[chimJunc]=$20;	
		gnNamesA[chimJunc]=$21;	
		gnNamesB[chimJunc]=$22;	
		gnTypesA[chimJunc]=$23;	
		gnTypesB[chimJunc]=$24;	
		juncSpanningReadIds[chimJunc]=$25;
		
		supportingPairsIds[chimJunc]="";
        inconsistentPairsIds[chimJunc]="";
		
        # Do all the possible gene pair ids combinations between the 2 list of genes 
        # overlapping junction donor and acceptor respectivelly.
        nbA=split(gnIdsA[chimJunc], gnListA, ",");
        nbB=split(gnIdsB[chimJunc], gnListB, ",");

        for (i=1;i<=nbA;i++)
        {
                for (j=1;j<=nbB;j++)
                {
                	if ((gnListA[i]!="") && (gnListB[j]!=""))
                	{
                		#print "eii", chimJunc, gnListA[i]"-"gnListB[j];
        			
        				donorGene=gnListA[i];
        				acceptorGene=gnListB[j];
        				
        				# Make for each DonorGene-AcceptorGene pair a list of chimeric junctions connecting them. 			
        				chimJunctionsList[donorGene"-"acceptorGene]=chimJunc","chimJunctionsList[gnListA[i]"-"gnListB[j]];  
        				nbSupportingPairs[donorGene"-"acceptorGene"::"chimJunc]=0;
						nbInconsistentPairs[donorGene"-"acceptorGene"::"chimJunc]=0;
        				supportingPairsIds[donorGene"-"acceptorGene"::"chimJunc]="";
        				inconsistentPairsIds[donorGene"-"acceptorGene"::"chimJunc]="";
        				   				
        				#print chimJunctionsList[gnListA[i]"-"gnListB[j]];
                	}
                }
        }		
        #print "##############";
	}	
}

{
	# print $0	

	## Get all relevant information about a discordant read pair and save into variables 
	pairId=$1;
	mapCoordMate1=$3;
	mapCoordMate2=$4;
	
	nbGnA=split($5, gnListA, ",");
    nbGnB=split($6, gnListB, ",");
    
    # print "numbers: ", nbGnA, nbGnB;
    
    # Do all the possible gene pair ids between the 2 list of genes 
    # overlapping discordant mate1 and mate2 respectivelly
    for (i=1;i<=nbGnA;i++)
    {
     	# print "number A: ", nbGnA;
     	# print "counter i: ", i;
    	for (j=1;j<=nbGnB;j++)
        {
        	# print "number B: ", nbGnB;
        	# print "counter j: ", j;
        	
        	# print gnListA[i]"-"gnListB[j];
        	
        	gnA=gnListA[i];
        	gnB=gnListB[j];
        	
        	## If both variables contain a gene id and the gene pair connection is supported by junction spanning reads 
        	if ((gnA!="") && (gnB != ""))
        	{ 
        		# A) If there are chimeric junctions connecting geneA and geneB (gene A is overlapping mate1 and gene B is overlapping mate2), respectively.
        		# Check if the discordant read-pair support them or not. 
        	
        		if (chimJunctionsList[gnA"-"gnB]!="")
        		{
        			donorGene=gnA;
        			acceptorGene=gnB;
        			
	        		mapCoordMateSupDonor=mapCoordMate1;
	        		mapCoordMateSupAcceptor=mapCoordMate2;	    	
	        		
	        		chimJunctions=chimJunctionsList[donorGene"-"acceptorGene];
	        		
	        	 	discordantPEsupport(pairId, donorGene, acceptorGene, mapCoordMateSupDonor, mapCoordMateSupAcceptor, chimJunctions);	
           		}
           		 
           		# B) If there are chimeric junctions connecting geneB and geneA (gene A is overlapping mate1 and gene B is overlapping mate2), respectively.
        		# Check if the discordant read-pair support them or not. 
           		if (chimJunctionsList[gnB"-"gnA] != "")
        		{
        			donorGene=gnB;
        			acceptorGene=gnA; 
        			
        			mapCoordMateSupDonor=mapCoordMate2;
	        		mapCoordMateSupAcceptor=mapCoordMate1;
	        		
	        		chimJunctions=chimJunctionsList[donorGene"-"acceptorGene];
	        		
	        	 	discordantPEsupport(pairId, donorGene, acceptorGene, mapCoordMateSupDonor, mapCoordMateSupAcceptor, chimJunctions);	
        		}
        	}        	
		}	
	}
}
END{
	for (chimJunc in chimJunctionsDic)
	{
		#print chimJunc, gnIdsA[chimJunc], gnIdsB[chimJunc];	
		
		maxSupportingPairs="0";
		maxInconsistentPairs="0";
		supPairsIds="na";
		inconPairsIds="na";
		
		# Do all the possible gene pair ids combinations between the 2 list of genes 
        # overlapping junction donor and acceptor respectivelly.
        nbA=split(gnIdsA[chimJunc], gnListA, ",");
        nbB=split(gnIdsB[chimJunc], gnListB, ",");

        for (i=1;i<=nbA;i++)
        {
        	for (j=1;j<=nbB;j++)
        	{
            	if ((gnListA[i]!="") && (gnListB[j]!=""))
               	{
             		donorGene=gnListA[i];		
             		acceptorGene=gnListB[j];
             		
                	#print "supporting: ", nbSupportingPairs[donorGene"-"acceptorGene"::"chimJunc];
                	#print "inconsistent: ", nbInconsistentPairs[donorGene"-"acceptorGene"::"chimJunc];
               	
               		if (nbSupportingPairs[donorGene"-"acceptorGene"::"chimJunc] > maxSupportingPairs)
               		{
               					maxSupportingPairs=nbSupportingPairs[donorGene"-"acceptorGene"::"chimJunc];
               					maxInconsistentPairs=nbInconsistentPairs[donorGene"-"acceptorGene"::"chimJunc];
               					supPairsIds=supportingPairsIds[donorGene"-"acceptorGene"::"chimJunc];
        						
        						if (inconsistentPairsIds[donorGene"-"acceptorGene"::"chimJunc] != "")
        						{
        							inconPairsIds=inconsistentPairsIds[donorGene"-"acceptorGene"::"chimJunc];			
               					}
               		}
               	}
            }
		}	
     	
     	## Compute total number of supporting paired-end (staggered+discordant) and
     	# percentage of inconsistent paired-end
     	totalNbPE=nbSpanningPE[chimJunc]+maxSupportingPairs;
		
		if ((maxInconsistentPairs == 0) && (maxSupportingPairs == 0))
		{
			percInconsistentPE="na"
		}	
		else	
		{
			percInconsistentPE=(maxInconsistentPairs/(maxInconsistentPairs+maxSupportingPairs))*100
		}
		
		## Print output
		row=chimJunc"\t"totalNbPE"\t"nbSpanningPE[chimJunc]"\t"nbStag[chimJunc]"\t"percStag[chimJunc]"\t"nbMulti[chimJunc]"\t"percMulti[chimJunc]"\t"maxSupportingPairs"\t"maxInconsistentPairs"\t"percInconsistentPE"\t"overlapA[chimJunc]"\t"overlapB[chimJunc]"\t"distExonBoundaryA[chimJunc]"\t"distExonBoundaryB[chimJunc]"\t"donorSS[chimJunc]"\t"acceptorSS[chimJunc]"\t"beg[chimJunc]"\t"end[chimJunc]"\t"sameChrStr[chimJunc]"\t"okGxOrder[chimJunc]"\t"dist[chimJunc]"\t"gnIdsA[chimJunc]"\t"gnIdsB[chimJunc]"\t"gnNamesA[chimJunc]"\t"gnNamesB[chimJunc]"\t"gnTypesA[chimJunc]"\t"gnTypesB[chimJunc]"\t"juncSpanningReadIds[chimJunc]"\t"supPairsIds"\t"inconPairsIds;
        
		print row;
	}
}

