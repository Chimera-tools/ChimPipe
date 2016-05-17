#!/usr/bin/env awk

# *****************************************************************************
	
#	bed12CorrectStrand.awk
	
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
# Takes a bed12 file and the variable "readDirectionality" as input and outputs a 6 fields bed file with the correct strand. If the variable readDirectionality is not defined it consider the data unstranded, so it does not change the strand. Example:
# awk -v readDirectionality=* -f bed12CorrectStrand.awk input.bed
# where * can be UNSTRANDED, MATE1_SENSE, MATE2_SENSE, SENSE OR ANTISENSE 

# example of input (".bed" text file)
##################

# chr1	14551	14601	SRR064286.5709675_HWI-EAS418:3:7:315:351_length=50#0/1	1	+	14551	14601	255,0,0	1	50	0
# chr1	14681	14731	SRR064286.5709675_HWI-EAS418:3:7:315:351_length=50#0/2	1	-	14681	14731	255,0,0	1	50	0
# chr1	14912	14962	SRR064286.7043546_HWI-EAS418:3:53:1606:39_length=50#0/1	1	+	14912	14962	255,0,0	1	50	0


## chr1    ChimPipe        alBlock 10525   10577   1       -       .       ReadName: "SRR201779.5012779_PATHBIO-SOLEXA2_30TUEAAXX:3:78:1387:1378_length=53#0/2";


BEGIN{OFS="\t"}
{
    chr=$1;
    beg=$2;
    end=$3;
    readId=$4;
    score=$5; # Number of hits	
    strand=$6;
    
    if (readDirectionality=="MATE1_SENSE")
    { 
		split(readId,id,"/"); 
		if(id[2]==2)
		{
			strand=(strand=="+" ? "-" : "+");
	 	} 
    }   	 
    else 	
    {	
		if (readDirectionality=="MATE2_SENSE")
		{
	   		split(readId,id,"/"); 
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
    
	print chr, beg, end, readId, score, strand;
}

     
     
