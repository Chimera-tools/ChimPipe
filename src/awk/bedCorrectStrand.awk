#!/usr/bin/env awk

# Takes a bed file and the variable "readDirectionality" as input and outputs a bed with the correct the strand. If the variable readDirectionality is not defined it consider the data unstranded, so it does not change the strand. Example:
# awk -v readDirectionality=* -f bedCorrectStrand.awk input.bed
# where * can be "MATE1_SENSE", "MATE2_SENSE", "MATE_STRAND_CSHL", "SENSE", "ANTISENSE" OR "NONE"

# example of input (".bed" text file)
##################
#chr1	3055369	3055470	HWI-ST985:73:C08BWACXX:8:2208:2017:40383/1	254	+	3055369	3055470	255,0,0	1	101	0
#chr1	3055453	3055554	HWI-ST985:73:C08BWACXX:8:2208:2017:40383/2	254	-	3055453	3055554	255,0,0	1	101	0
#chr1	3068681	3068782	HWI-ST985:73:C08BWACXX:8:1206:4407:71813/1	180	+	3068681	3068782	255,0,0	1	101	0
#chr1	3068708	3068809	HWI-ST985:73:C08BWACXX:8:1206:4407:71813/2	180	-	3068708	3068809	255,0,0	1	101	0

BEGIN{OFS="\t"}
{
    if (readDirectionality=="MATE1_SENSE")
    { 
	split($4,id,"/"); 
	if(id[2]==2)
	{
	    $6=($6=="+" ? "-" : "+");
	} 
    }   	 
    else 	
    {	
	if (readDirectionality=="MATE2_SENSE")
	{
	    split($4,id,"/"); 
	    if(id[2]==1)
	    {
		    $6=($6=="+" ? "-" : "+");
	    } 
	}
	else 
	{
	    if (readDirectionality=="ANTISENSE")
	    {
		$6=($6=="+" ? "-" : "+"); 
	    }
	}
    }
    print $0;	
}

     
     
