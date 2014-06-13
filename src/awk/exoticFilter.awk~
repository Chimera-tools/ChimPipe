#!/usr/bin/env awk

# Takes as input a split mapping gem file from gem-rna-mapper (June 2013) and outputs reads mapping uniquely, in 2 blocks and in different chromosome or strand. 

# example of input
##################
# SINATRA_0006:1:1:7430:930#0/1	NCCTTCTCTTCGCTCCTGGTGTAAGGTATGGTACATAAGAGTCCAATGCTATTTGCGCAAGTGCTAGGGTAACGAG	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB0:0:0:0:1	chr1:+:80324753:15::chr12:+:112677685:61
# SINATRA_0006:1:1:3790:968#0/1	NAACTCATCATAGTGTTCCTGCATCTCCACATCGCTCACGGCACAGTGTGAGCCGTCAGCCGTCTGTGCACTGTTT	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB0:0:1	chr21:+:44515810:44::chr21:+:44521476:32



$NF!="-"{
    n=split($NF,a,","); 
    if(n==1)
    {
		n2=split($NF,b,"::"); 
		if(n2==2)
		{
	   		split(b[1],b1,":"); 
	    	split(b[2],b2,":"); 
	    	if (stranded==""||stranded==0)
	    	{	
	    		if ((b1[1]!=b2[1]) || (b1[2]!=b2[2]))
	    		{
	    			print $0
	    		}
			}
			else
			{
				if ((b1[1]!=b2[1]) || (b1[2]!=b2[2]) || ((b1[2]=="+") && (b2[2]=="+") && (b1[3]>b2[3])) || ((b1[2]=="-") && (b2[2]=="-") && (b2[3]>b1[3])))	    
				{
		    		print $0;
				}
			}
		}
    }
}

