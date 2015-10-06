#!/usr/bin/env awk

# *****************************************************************************
	
#	gffCorrectStrand.awk
	
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
# Takes a gff file and the variable "readDirectionality" as input and outputs a gff with the correct the strand. If the variable readDirectionality is not defined it consider the data unstranded, so it does not change the strand. Example:
# awk -v readDirectionality=* -f bedCorrectStrand.awk input.bed
# where * can be UNSTRANDED, MATE1_SENSE, MATE2_SENSE, SENSE OR ANTISENSE 

# example of input (".gff" text file)
##################

## chr1    ChimPipe        alBlock 10525   10577   1       -       .       ReadName: "SRR201779.5012779_PATHBIO-SOLEXA2_30TUEAAXX:3:78:1387:1378_length=53#0/2";



BEGIN{OFS="\t"}
{
    strand=$7;
    split($10,read,"\"");
    
    if (readDirectionality=="MATE1_SENSE")
    { 
		split(read[2],id,"/"); 
		if(id[2]==2)
		{
			strand=(strand=="+" ? "-" : "+");
	 	} 
    }   	 
    else 	
    {	
		if (readDirectionality=="MATE2_SENSE")
		{
	   		split(read[2],id,"/"); 
	    	if(id[2]==1)
	    	{
		    	strand=(strand=="+" ? "-" : "+");
	    	} 
		}
		else 
		{
	    	if (readDirectionality=="ANTISENSE")
	    	{
				strand=(strand=="+" ? "-" : "+");
	    	}
		}
    }
    
	print $1, $2, $3, $4, $5, $6, strand, $8, $9, $10;
}

     
     
