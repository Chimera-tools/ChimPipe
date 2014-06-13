

# chimjunc_to_gff.awk

# junc=/users/rg/brodriguez/Chimeras_project/Chimeras_detection_pipeline/Chimsplice/Versions/V0.4.0/Test/Intermediate_files/distinct_junctions_nbstaggered_nbtotalsplimappings_minBeg_maxEnd_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn_LID16627.txt 
# chr7_128502985_+:chr1_231556974_- 1 4 128502926 231556959
# 10777 (5 fields)   *** juncid, nbstag, nbtot, beg, end (from staggered_to_junction.awk)

{
    split($1,a,":"); 
    split(a[1],a1,"_"); 
    split(a[2],a2,"_"); 
    
    if (stranded==1) # stranded data 
    {	
    	if (a1[2]<$4)	# part 1
    	{
			print a1[1], ".", "exon", a1[2], $4, ".", a1[3], ".", "junc:", $1; 
    	}
    	else
    	{
			print a1[1], ".", "exon", $4, a1[2], ".", a1[3], ".", "junc:", $1;     
    	}
    
    	if (a2[2]<$5)	# part 2
    	{
			print a2[1], ".", "exon", a2[2], $5, ".", a2[3], ".", "junc:", $1; 
    	}
    	else
    	{
			print a2[1], ".", "exon", $5, a2[2], ".", a2[3], ".", "junc:", $1;     
    	}
    }
    
    else  # unstranded data 
    {	
    	if (a1[2]<$4)	# part 1
    	{
			print a1[1], ".", "exon", a1[2], $4, ".", ".", ".", "junc:", $1; 
    	}
    	else
    	{
			print a1[1], ".", "exon", $4, a1[2], ".", ".", ".", "junc:", $1;     
    	}
    
    	if (a2[2]<$5)	# part 2
    	{
			print a2[1], ".", "exon", a2[2], $5, ".", ".", ".", "junc:", $1; 
    	}
    	else
    	{
			print a2[1], ".", "exon", $5, a2[2], ".", ".", ".", "junc:", $1;     
    	}
    }
}    
    
#    if (stranded==0) 
#    {
#		print a1[1], ".", "exon", $4, a1[2], ".", ".", ".", "junc:", $1; 
#		print a2[1], ".", "exon", a2[2], $5, ".", ".", ".", "junc:", $1;
#    } 
#    else 
#    {
#		if ((a1[3]=="+") && (a2[3]=="-"))
#		{
#		    print a1[1], ".", "exon", $4, a1[2], ".", a1[3], ".", "junc:", $1; 
#		    print a2[1], ".", "exon", $5, a2[2], ".", a2[3], ".", "junc:", $1;
#		} 
#		else 
#		{
#		    if ((a1[3]=="-") && (a2[3]=="+"))
#		    {
#				print a1[1], ".", "exon", a1[2], $4, ".", a1[3], ".", "junc:", $1; 
#				print a2[1], ".", "exon", a2[2], $5, ".", a2[3], ".", "junc:", $1;
#		    }
#		    else
#		    {
#				print a1[1], ".", "exon", $4, a1[2], ".", a1[3], ".", "junc:", $1; 
#				print a2[1], ".", ".", a2[2], $5, ".", a2[3], ".", "junc:", $1;
#		    }
#		}
#    }

