#!/usr/bin/env awk

# Takes a split mapping gem file from the gem-rna-mapper (June 2013) and the variable "readDirectionality" as input and outputs the unique 2 block split mappings with the correct strand. If the variable readDirectionality is not defined it consider the data unstranded, so it does not change the strand. Example:
# awk -v readDirectionality=* -f bedCorrectStrand.awk input.map
# where * can be "MATE1_SENSE", "MATE2_SENSE", "MATE_STRAND_CSHL", "SENSE", "ANTISENSE" OR "NONE"

# example of input (".map" file)
##################
# SINATRA_0006:1:1:7430:930#0/1	NCCTTCTCTTCGCTCCTGGTGTAAGGTATGGTACATAAGAGTCCAATGCTATTTGCGCAAGTGCTAGGGTAACGAG	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB 0:0:0:0:1	chr1:+:80324753:15::chr12:+:112677685:61
# SINATRA_0006:1:1:3790:968#0/1	NAACTCATCATAGTGTTCCTGCATCTCCACATCGCTCACGGCACAGTGTGAGCCGTCAGCCGTCTGTGCACTGTTT	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB 0:0:1	chr21:+:44515810:44::chr21:+:44521476:32

BEGIN{OFS="\t"}
$NF!="-"{
	rev=0;
    n=split($NF,a,","); 
    if(n==1)
    {
		n2=split($NF,b,"::"); 
		if(n2==2)
		{
	    	split(b[1],b1,":"); 
	    	split(b[2],b2,":");
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
	    		print $1, $2, $3, $4, b1[1]":"strA":"b1[3]":"b1[4]"::"b2[1]":"strB":"b2[3]":"b2[4];
			}
			else
			{
				print $1, $2, $3, $4, b2[1]":"strB":"b2[3]":"b2[4]"::"b1[1]":"strA":"b1[3]":"b1[4];
			}
		}
    }
}	

