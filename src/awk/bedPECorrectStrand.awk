#!/usr/bin/env awk

# Takes a bed file and the variable "readDirectionality" as input and outputs a bed with the correct the strand. If the variable readDirectionality is not defined it consider the data unstranded, so it does not change the strand. Example:
# awk -v readDirectionality=* -f bedCorrectStrand.awk input.bed
# where * can be "MATE1_SENSE", "MATE2_SENSE", "MATE_STRAND_CSHL", "SENSE", "ANTISENSE" OR "NONE"

# example of input (".bed" text file)
##################

## chr1 752991 753018 chr1 754245 754250 HWI-ST985:73:C08BWACXX:8:2208:2017:40383/1  1000 + +



BEGIN{OFS="\t"}
{
    strA=$9;
	strB=$10;
	rev=0;
    if (readDirectionality=="MATE1_SENSE")
    { 
		split($7,id,"/"); 
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
	   		split($7,id,"/"); 
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
		print $1, $2, $3, $4, $5, $6, $7, $8, strA, strB;
	}
	else
	{
		print $4, $5, $6, $1, $2, $3, $7, $8, strB, strA;
	}
}

     
     
