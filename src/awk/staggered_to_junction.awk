
# staggered_to_junction.awk

# example
# awk -v stranded=$stranded $STAG $outdir/distinct_staggered_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt 

# input
# chr10_101153356_101153370_+:chr16_4660515_4660552_- 2
# chr10_102122971_102123002_-:chr14_75370133_75370153_+ 4

# output

{
    split($1,a,":"); 
    split(a[1],a1,"_"); 
    split(a[2],a2,"_"); 
    if (stranded==1)  # stranded data
    {
		if ((a1[4]==a2[4]))  # same strand
		{
			if ((a1[4]=="+")) # part1 + and part2 +
		    {
		    	junc=a1[1]"_"a1[3]"_"a1[4]":"a2[1]"_"a2[2]"_"a2[4];
		    	nbStag[junc]++;
		    	nbTotal[junc]=nbTotal[junc]+$2;
		    	if ((beg[junc]=="") || (beg[junc]>a1[2]))
		    	{
		    		beg[junc]=a1[2];	
		    	}
		    	if ((end[junc]=="") || (end[junc]<a2[3]))
		    	{
		    		end[junc]=a2[3];
		    	}
		    }
		    else # part1 - and part2 - (Reverse the blocks to write them in biological order)
		    {
		    	junc=a2[1]"_"a2[2]"_"a2[4]":"a1[1]"_"a1[3]"_"a1[4];
		    	nbStag[junc]++;
		    	nbTotal[junc]=nbTotal[junc]+$2;
		    	if ((beg[junc]=="") || (beg[junc]<a2[3]))
		    	{
		    		beg[junc]=a2[3];
		    	}		    	
		    	if ((end[junc]=="") || (end[junc]>a1[2]))
		    	{
		    		end[junc]=a1[2];	
		    	}
		    	
		    }	
		}
		else	# diff strand
		{
			if ((a1[4]=="+") && (a2[4]=="-")) # part1 + and part2 -
			{
			    junc=a1[1]"_"a1[3]"_"a1[4]":"a2[1]"_"a2[3]"_"a2[4];
			    nbStag[junc]++;
			    nbTotal[junc]=nbTotal[junc]+$2;
			    if ((beg[junc]=="") || (beg[junc]>a1[2]))
		    	{
		    		beg[junc]=a1[2];	
		    	}
		    	if ((end[junc]=="") || (end[junc]>a2[2]))
		    	{
		    		end[junc]=a2[2];
		    	}
			}
			else	# part1 - and part2 +
			{
			    if ((a1[4]=="-") && (a2[4]=="+")) # Part1 - and part2 +
			   	{
					junc=a1[1]"_"a1[2]"_"a1[4]":"a2[1]"_"a2[2]"_"a2[4];
			    	nbStag[junc]++;
			    	nbTotal[junc]=nbTotal[junc]+$2;
			    	if ((beg[junc]=="") || (beg[junc]<a1[3]))
		    		{
		    			beg[junc]=a1[3];	
		    		}
		    		if ((end[junc]=="") || (end[junc]<a2[3]))
		    		{
		    			end[junc]=a2[3];
		    		}
			   	}
			}
		}		
    }
    else	# Unstranded data
    {
		if (a1[4]==a2[4])	# same strand 
		{
	    	if (a1[1]==a2[1])	# same chr  
	    	{
				junc=a1[1]"_"a1[3]":"a2[1]"_"a2[2];
				nbStag[junc]++;
				nbTotal[junc]=nbTotal[junc]+$2;
				if ((beg[junc]=="") || (beg[junc]>a1[2]))
		    	{
		    		beg[junc]=a1[2];	
		    	}
		    	if ((end[junc]=="") || (end[junc]<a2[3]))
		   	 	{
		    		end[junc]=a2[3];
		    	}
	    	}
	    	else    # diff chromosome
	    	{
	    		if (a1[4]=="+")	# part1 + and part2 + 
	    		{	    	
					if (a1[2]<=a2[2]) # numerical order 	
					{
		    			junc=a1[1]"_"a1[3]":"a2[1]"_"a2[2];
						nbStag[junc]++;
						nbTotal[junc]=nbTotal[junc]+$2;
						if ((beg[junc]=="") || (beg[junc]>a1[2]))
		    			{
		    				beg[junc]=a1[2];	
		    			}
		    			if ((end[junc]=="") || (end[junc]<a2[3]))
		    			{
		    				end[junc]=a2[3];
		    			}
					}
					else	# not numerical order -> write the junctions in numerical order
					{
		    			junc=a2[1]"_"a2[2]":"a1[1]"_"a1[3];
						nbStag[junc]++;
						nbTotal[junc]=nbTotal[junc]+$2;
						if ((beg[junc]=="") || (beg[junc]<a2[3]))
		    			{
		    				beg[junc]=a2[3];	
		    			}
		    			if ((end[junc]=="") || (end[junc]>a1[2]))
		    			{
		    				end[junc]=a1[2];
		    			}			
					}	
		    	}
		    	else	# part1 - and part2 -
		    	{
		    		if (a1[2]<=a2[2]) # numerical order 
					{
		    			junc=a1[1]"_"a1[2]":"a2[1]"_"a2[3];
						nbStag[junc]++;
						nbTotal[junc]=nbTotal[junc]+$2;
						if ((beg[junc]=="") || (beg[junc]<a1[3]))
		    			{
		    				beg[junc]=a1[3];	
		    			}
		    			if ((end[junc]=="") || (end[junc]>a2[2]))
		    			{
		    				end[junc]=a2[2];
		    			}
					}
					else	# not numerical order -> write the junctions in numerical order
					{	
		    			junc=a2[1]"_"a2[3]":"a1[1]"_"a1[2];
						nbStag[junc]++;
						nbTotal[junc]=nbTotal[junc]+$2;
						if ((beg[junc]=="") || (beg[junc]>a2[2]))
		    			{
		    				beg[junc]=a2[2];	
		    			}
		    			if ((end[junc]=="") || (end[junc]<a1[3]))
		    			{
		    				end[junc]=a1[3];
		    			}			
					}	
		    	}					
	    	}
		}
		else 	# diff strand 
		{
	    	if(a1[4]=="+")  # part1 + and part2 - 
	    	{
				if (a1[2]<=a2[2]) # numerical order 
				{
		    		junc=a1[1]"_"a1[3]":"a2[1]"_"a2[3];
					nbStag[junc]++;
					nbTotal[junc]=nbTotal[junc]+$2;
					if ((beg[junc]=="") || (beg[junc]>a1[2]))
		    		{
		    			beg[junc]=a1[2];	
		    		}
		    		if ((end[junc]=="") || (end[junc]>a2[2]))
		    		{
		    			end[junc]=a2[2];
		    		}
				}
				else	# not numerical order -> write the junctions in numerical order
				{
		    		junc=a2[1]"_"a2[3]":"a1[1]"_"a1[3];
					nbStag[junc]++;
					nbTotal[junc]=nbTotal[junc]+$2;
					if ((beg[junc]=="") || (beg[junc]>a2[2]))
		    		{
		    			beg[junc]=a2[2];	
		    		}
		    		if ((end[junc]=="") || (end[junc]>a1[2]))
		    		{
		    			end[junc]=a1[2];
		    		}			
				}
	    	}
	    	else  # part1 - and part2 +
	    	{
				if(a1[2]<=a2[2])	# numerical order
				{
		    		junc=a1[1]"_"a1[2]":"a2[1]"_"a2[2];
		    		nbStag[junc]++;
		    		nbTotal[junc]=nbTotal[junc]+$2;
		    		if ((beg[junc]=="") || (beg[junc]<a1[3]))
		    		{
		    			beg[junc]=a1[3];	
		    		}
		    		if ((end[junc]=="") || (end[junc]<a2[3]))
		    		{
		    			end[junc]=a2[3];
		    		}
				}
				else	# not numerical order -> write the junctions in numerical order
				{
		    		junc=a2[1]"_"a2[2]":"a1[1]"_"a1[2];
		    		nbStag[junc]++;
		    		nbTotal[junc]=nbTotal[junc]+$2;
		    		if ((beg[junc]=="") || (beg[junc]<a2[3]))
		    		{
		    			beg[junc]=a2[3];
		    		}
		    		if ((end[junc]=="") || (end[junc]<a1[3]))
		    		{
		    			end[junc]=a1[3];
		    		}
				}
	    	}
		}
    }
} 
END {
	for (junc in nbStag)
	{
		print junc, nbStag[junc], nbTotal[junc], beg[junc], end[junc];
	}
}


