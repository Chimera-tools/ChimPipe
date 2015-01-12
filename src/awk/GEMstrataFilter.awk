#!/usr/bin/env awk

# *****************************************************************************
	
#	GEMstrataFilter.awk
	
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
# ******************************************************************************

# Description
###############
# Takes as input a file with alignments in GEM format and select only the alignments from the best tratum. 

##  Input (GEM file)
#SRR064438.11963715_HWI-EAS418:6:79:787:1144_length=50#0 GCCAGGCTAGGAGGAATACCTGTGGGAGTTGTTGCTGTAGAAACCCGAAC GCTGGAGCTTGGCTTCAGAATCCAGGTTTGCTGGATCAGCTGGGATACTT   BCCCCBCCBCBBCBBCCBBBBBABBB@BABBBBBBBB@ABABA?BB@:?	< BCCCCBBCBBBBBBBBBBBBBBBABB@BBBABBBABBBBBBBBB?AA=B@    0:0:0:0:1:2		chr17:-:35479505:5>3069*45::chr17:+:35479447:2>1+1>2-1C43:::16000,chr17:-:35479505:AT48::chr17:+:35478403:6>3+T>1040*43:::12032,chr17:-:35479505:AT48::chr17:+:35479447:2>1+1>2-1C43:::12160

##  Output (filtered GEM file with only one strata)
#SRR064438.11963715_HWI-EAS418:6:79:787:1144_length=50#0 GCCAGGCTAGGAGGAATACCTGTGGGAGTTGTTGCTGTAGAAACCCGAAC GCTGGAGCTTGGCTTCAGAATCCAGGTTTGCTGGATCAGCTGGGATACTT   BCCCCBCCBCBBCBBCCBBBBBABBB@BABBBBBBBB@ABABA?BB@:?	< BCCCCBBCBBBBBBBBBBBBBBBABB@BBBABBBABBBBBBBBB?AA=B@    0:0:0:0:2     chr17:-:35479505:5>3069*45::chr17:+:35479447:2>1+1>2-1C43:::16000


# Usage:
########
# awk -v OFS="\t" <INPUT_ARGS> -f SAMfilter.awk mappings.map



# Print if unmapped
$NF=="-"{
	print $0;
}

# If mapped
$NF!="-"{
	stratum=0; # count 
	matchSummaryOut=""; 
	alignmentsIn=$NF; 
	
	# If NOT paired alignment
	if (NF==5)
	{
		matchSummaryIn=$4;
	}
	else # If paired alignment
	{
		matchSummaryIn=$6
	}
	 
	gsub(/+/, ":", matchSummaryIn);  # In some cases the match summary contains "+" instead of ":". Replace the "+" by ":" to have all the stratum separated with ":". Ejem: 0:0:1
	
	maxNbStrata=split(matchSummaryIn,strata,":"); 
	maxNbAlignments=split(alignmentsIn,alignments,","); 
	
	# Parse the match summary from first to last stratum to take only the best statum
	while (stratum++<maxNbStrata) 
	{
		if (strata[stratum]=="0")
		{
			matchSummaryOut="0:"matchSummaryOut
		}
		else
		{
			maxNbReportedAlignments=strata[stratum];
			matchSummaryOut=(matchSummaryOut)(strata[stratum]); 
			break;
		}
	} 
		
	if (maxNbReportedAlignments==maxNbAlignments)
	{
		alignmentsOut=$NF;
	}			
	else
	{
		alignment=0;
		alignmentsOut=""; 
		while (alignment++<maxNbReportedAlignments)
		{
			if (alignmentsOut == "")
			{
				alignmentsOut=(alignments[alignment]);
			}
			else
			{
				alignmentsOut=(alignmentsOut)","(alignments[alignment]);		
			}
		}
	}
	
	# Replace the previous alignment and match summary with the new one. 
	if (NF==5)
	{
		matchSummaryIn=$4;
	}
	else # If NOT paired alignment
	{
		matchSummaryIn=$6
	}
	
	$NF=alignmentsOut
	$(NF-1)=matchSummaryOut
	
	# PRINT OUT THE FILTERED ALIGNMENTS
	# If NOT paired alignment
	if (NF==5)
	{
		print $1, $2, $3, $4, $5; 
		
	}
	else # If paired alignment
	{
		print $1, $2" "$3, $4" "$5, $6, $7;
	}
}



