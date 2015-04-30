#!/usr/bin/env awk

# *****************************************************************************
	
#	add_mateInfo_SAM.awk
	
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

# Takes as input a SAM file with paired-end mapped reads, tags the read ids for mate 1 and mate 2 as /1 and /2, respectively, and print the SAM to the stdout. 

##  Input (file in SAM alignment format)

# SRR201779.2057082_PATHBIO-SOLEXA2_30TUEAAXX:3:31:1009:1024_length=53#0	163	chr1	877861	254	8M70N45M	=	878079	271	CTGGCCCGGCTGGAGCTGCCCGCCGACCTCCTGCGGCAGAAGGAGCTGGAGAG	CBBBCBCBBBBBB;C@;?BCCBBBBBBBBBBBBBBBBBBBBBBBB@;A>@ABB	
# SRR201779.2057082_PATHBIO-SOLEXA2_30TUEAAXX:3:31:1009:1024_length=53#0	163	chr1	878392	254	53M	=	878753	734	ACTGCCCCTGGGCTTCCCTTATGCCGTCAGCCCCTACTTCCACACAGGCGCGG	BBBBBBCBBBBBBBBBCCBBBBBCBBBBBBBCCB@;ABBB@6?2?2?BBBBB9	

##  Output (file in SAM alignment format with the read ids tagged with /1 and /2 for mate 1 and mate 2 respectivelly)

# SRR201779.2057082_PATHBIO-SOLEXA2_30TUEAAXX:3:31:1009:1024_length=53#0/1	163	chr1	877861	254	8M70N45M	=	878079	271	CTGGCCCGGCTGGAGCTGCCCGCCGACCTCCTGCGGCAGAAGGAGCTGGAGAG	CBBBCBCBBBBBB;C@;?BCCBBBBBBBBBBBBBBBBBBBBBBBB@;A>@ABB	
# SRR201779.2057082_PATHBIO-SOLEXA2_30TUEAAXX:3:31:1009:1024_length=53#0/2	163	chr1	878392	254	53M	=	878753	734	ACTGCCCCTGGGCTTCCCTTATGCCGTCAGCCCCTACTTCCACACAGGCGCGG	BBBBBBCBBBBBBBBBCCBBBBBCBBBBBBBCCB@;ABBB@6?2?2?BBBBB9	

# Usage example:
################
# awk -v OFS="\t" -f add_mateInfo_sam.awk mappings.sam

{
	if (($1 ~ /^@/)||($1 ~ /\/1$/)||($1 ~ /\/2$/))  # Skip the header and read ids already tagged 
	{	
		print $0;
	}
	else
	{
		if (and($2,0x40)) # If the SAM Bitwise flag specifies "first segment in the template" tag the read id with /1
		{
			$1=$1"\/1";
		}
		else if (and($2,0x80)) # If SAM Bitwise specifies "last segment in the template" tag the read id with /2
		{
			$1=$1"\/2";
		} 
		print $0;
	}
}

