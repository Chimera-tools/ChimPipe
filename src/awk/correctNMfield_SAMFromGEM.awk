#!/usr/bin/env awk

# *****************************************************************************
	
#	correctNMfield_SAMfromGEM.awk
	
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

# Takes as input a SAM file. Check if it was produced from a GEM mapping file with gem-2-sam and if it is the case corrects the edit distance to the reference in the NM:i:X optional field. Needed because gem-2-sam considers the skipped regions from the reference (N in the CIGAR string) as mismatches. I think it should not be included since it adds mismatches to reads spanning exon-exon junctions, while this is not really a mismatch and another mappers like tophat does not consider them as mismatches.  

##  Input (file in SAM format)

# SRR201779.5153708_PATHBIO-SOLEXA2_30TUEAAXX:3:80:590:450_length=53#0	99	chr1	883572	179	41M257N12M	=	886508	2989	TACTCCAGGGTGAGGTCGTACAGCTGCTCCACCACGCCGCCCCAGAACGCCTT	BBBBB?;<>9039;BB@;AA<99927<9A9&1;9&,?<63673&0.36&3?99	RG:Z:0	NH:i:1	NM:i:3	XT:A:U	md:Z:34G4T1>257*2G1T7

# SRR201779.1181413_PATHBIO-SOLEXA2_30TUEAAXX:3:18:686:365_length=53#0     355     chr1    957839  125     4M12814N48M5340N1M      =       976158  18372   *       *       RG:Z:0  NH:i:4  NM:i:2  XT:A:U	md:Z:4>12814*48>5340*1

##  Output (file in SAM format with the correct number of mismatches in the NM:i:X optional field)

# SRR201779.5153708_PATHBIO-SOLEXA2_30TUEAAXX:3:80:590:450_length=53#0	99	chr1	883572	179	41M257N12M	=	886508	2989	TACTCCAGGGTGAGGTCGTACAGCTGCTCCACCACGCCGCCCCAGAACGCCTT	BBBBB?;<>9039;BB@;AA<99927<9A9&1;9&,?<63673&0.36&3?99	RG:Z:0	NH:i:1	NM:i:2	XT:A:U	md:Z:34G4T1>257*2G1T7

# SRR201779.1181413_PATHBIO-SOLEXA2_30TUEAAXX:3:18:686:365_length=53#0     355     chr1    957839  125     4M12814N48M5340N1M      =       976158  18372   *       *       RG:Z:0  NH:i:4  NM:i:2  XT:A:U	md:Z:4>12814*48>5340*1


# Usage example:
################
# awk -v OFS="\t" -f correctNMfield_SAMfromGEM.awk mappings.sam

 
/^@/{ 	     # Print the header
	if ($0 ~ "PG:GEM") # Check if reads mapped with GEM
	{
		GEM="1";
	}
	print $0;
}

! /^@/{   # Print alignments 
	
	if ((GEM=="1")&&( ! and($2,0x4)))  # If read split-mapped with GEM correct the number of mismatches field 
	{
		n=split($6,cigar,"N");  # Count the number of Ns in the CIGAR
		nbSkipped=n-1; 
		split($14,NM,":");
		$14="NM:i:"(NM[3]-nbSkipped);  # Substract the nb of Ns from the nb of mismatches to compute the real nb. of mismatches
	}
	print $0;
}

