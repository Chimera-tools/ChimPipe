#!/usr/bin/env awk

# *****************************************************************************
	
#	gemCorrectStrand.awk
	
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

# Takes a split mapping gem file from the gem-rna-mapper (June 2013) and the variable "readDirectionality" as input and outputs the unique 2 block split mappings with the correct strand and block orientation. If the variable readDirectionality is not defined it consider the data unstranded, so it does not change the strand. Example:
# awk -v readDirectionality=* -f bedCorrectStrand.awk input.map
# where * can be UNSTRANDED, MATE1_SENSE, MATE2_SENSE, SENSE OR ANTISENSE 

# example of input (".map" file)
##################
# SINATRA_0006:1:1:7430:930#0/1	NCCTTCTCTTCGCTCCTGGTGTAAGGTATGGTACATAAGAGTCCAATGCTATTTGCGCAAGTGCTAGGGTAACGAG	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB 0:0:0:0:1	chr1:+:80324753:15::chr12:+:112677685:61
# SINATRA_0006:1:1:3790:968#0/1	NAACTCATCATAGTGTTCCTGCATCTCCACATCGCTCACGGCACAGTGTGAGCCGTCAGCCGTCTGTGCACTGTTT	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB 0:0:1	chr21:+:44515810:44::chr21:+:44521476:32

# Ouput if MATE2_SENSE data
###########################
# SINATRA_0006:1:1:7430:930#0/1	NCCTTCTCTTCGCTCCTGGTGTAAGGTATGGTACATAAGAGTCCAATGCTATTTGCGCAAGTGCTAGGGTAACGAG	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB 0:0:0:0:1	chr12:-:112677685:61::chr1:-:80324753:15
# SINATRA_0006:1:1:3790:968#0/1	NAACTCATCATAGTGTTCCTGCATCTCCACATCGCTCACGGCACAGTGTGAGCCGTCAGCCGTCTGTGCACTGTTT	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB 0:0:1	chr21:-:44521476:32::chr21:-:44515810:44

# Usage example:
# awk -v readDirectionality="MATE2_SENSE" -f gemCorrectStrand.awk

BEGIN{OFS="\t"}
$NF!="-"{
	rev=0;
    mappings="";
    
    #n=split($NF,a,","); 
    nbHits=split($NF,alignments,","); 
    for (align in alignments)
    {
		nbBlocks=split(alignments[align],blocks,"::"); 
		if(nbBlocks==2)
		{
	   	 	split(blocks[1],b1,":"); 
	    	split(blocks[2],b2,":");
	    	strA=b1[2];
	    	strB=b2[2];
	    	if (readDirectionality=="MATE1_SENSE")
	    	{ 
				split($1,id,"/"); 
				if(id[2]==2)
				{	
			    	strA=(strA=="+" ? "-" : "+");
			    	strB=(strB=="+" ? "-" : "+");
			    	rev="1";		   		    
				} 
	    	}   	 
	    	else 
	    	{
				if (readDirectionality=="MATE2_SENSE")
				{
				    split($1,id,"/"); 
				    if(id[2]==1)
				    {
						strA=(strA=="+" ? "-" : "+");
						strB=(strB=="+" ? "-" : "+");
						rev="1";	   				
				    } 
				}
    			else 
				{
				    if (readDirectionality=="ANTISENSE")
				    {
						strA=(strA=="+" ? "-" : "+");
						strB=(strB=="+" ? "-" : "+");
						rev="1";
				    }
				}
	    	}
	    	if (rev=="0")
	    	{	
	    		 map=b1[1]":"strA":"b1[3]":"b1[4]"::"b2[1]":"strB":"b2[3]":"b2[4];
			}
			else
			{
				 map=b2[1]":"strB":"b2[3]":"b2[4]"::"b1[1]":"strA":"b1[3]":"b1[4];
			}
			
			if (mappings == "")
			{
				mappings=map
			}
			else
			{
				mappings=mappings","map
			}
		}
    }
    print $1, $2, $3, $4, mappings;
}	

