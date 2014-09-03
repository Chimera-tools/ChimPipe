#!/usr/bin/env awk

# Takes as input a file containing staggered reads split-mapping in two exons from different genes and the 8-kmer on each side of the blocks. It uses this information to make the correct chimeric junctions joining the donnor and acceptor sides of the blocks, compute the number of staggered reads supporting the junctions, their maximum beginning and end coordinates and their consensus splice sites. 

# input
# chr10_100013508_100013530_-:chr3_50392315_50392367_-	1	CTACTTGG	TGATGAAC	GTTTCCCT	GAGCCGAC		
# chr10_100177011_100177032_-:chr2_85133614_85133667_+	3	CTCTGCTG	TTATCCTG	AGAGCCAC	GATGCCAC		
# chr10_101464231_101464286_+:chr21_41447039_41447058_-	1	AAGGTGTC	CTTTGCCT	GGGTCTTG	CTGGTGAC		 
# chr10_101874926_101874942_-:chr17_38137273_38137331_+	4	CTGCGCGC	CCGCAGCA	GGCGAAAC	ATGAAAGA		 
# chr10_102045932_102045960_-:chrX_63946916_63946962_-	2	CTGCTGAA	GATGCCGC	ATGGACAT	AAAGGAAC		 

# output
# chr3_50392367_+:chr10_100013508_+ 1 1 50392315 100013530 GT AG
# chr2_85133614_-:chr10_100177011_+ 2 6 85133671 100177032 GT AG
# chr21_41447058_+:chr10_101464286_- 1 1 41447039 101464231 GT AG
# chr17_38137273_-:chr10_101874926_+ 1 4 38137331 101874942 GT AG
# chrX_63946962_+:chr10_102045932_+ 1 2 63946916 102045960 GT AG

# Usage example:
# awk -v strandedness="1" -v CSS='GT+AG,GC+AG,ATATC+A.,GTATC+AT' -f $MAKE_JUNCTIONS $outdir/staggered_nb_dnt.txt


# How the script find the position of the splice sites regarding on the strand the blocks map?? 
###############################################################################################
# Their position depend on the block orientation and this depend on last term on the strand the read maps.
##########################################################################################################
# +/+ 
#	  BLOCK1			BLOCK2
# c1a-------c1b		c2a-------c2b 
# The junction should be made between the internal coordinates c1b and c2a. 
 
# c1a-------c1b	ss1 ss2 c2a-------c2b 

# c1a-------c1b	ss2 ss1 c2a-------c2b				

# +/-
#	  BLOCK1			BLOCK2
# c1a-------c1b		c2a-------c2b 
# The junction should be made between the coordinates c1b and c2b. .
 
# c1a-------c1b	ss1  c2a-------c2b ss2

# c1a-------c1b	ss2  c2a-------c2b ss1				

# -/+ 
#	  BLOCK1			BLOCK2
# c1a-------c1b		c2a-------c2b 
# The junction should be made between the coordinates c1b and c2b. 
 
# ss1 c1a-------c1b	  ss2 c2a-------c2b 

# ss2 c1a-------c1b	  ss1 c2a-------c2b 				

# -/- 
#	  BLOCK1			BLOCK2
# c1a-------c1b		c2a-------c2b 
# The junction should be made between the coordinates c1a and c2b, since the blocks are in bam convention (genomic order instead biological one, so for -- the order of the blocks is the opposite than for ++). 
 
# ss1 c1a-------c1b	c2a-------c2b ss2

# ss2 c1a-------c1b	c2a-------c2b ss1		



function revComp(seq)
{
     s="";
     n=split(seq, a, "");
     for (x=n; x >= 1; x--)
     {
     	
     	if (a[x]=="A")
     	{
     		comp="T";
     	}
     	else if (a[x]=="T")
     	{
     		comp="A";
     	}
     	else if (a[x]=="G")
     	{
     		comp="C";
     	}
     	else if (a[x]=="C")
     	{
     		comp="G";
     	}
     	else
     	{
     		comp="."
     	}
     	s=(s)(comp);
     }
     return s
}

function findStrand(seq1,seq2,mapStr1,mapStr2,CSS)
{
	# Initializing variables
	juncStr1="";
	juncStr2="";
	juncDonor="";
	juncAcceptor="";
	
	nbcss=split(CSS,b,",");
	for (i = 1; i <= nbcss; i++) # Iterate over the consensus splice sites given as input
	{
		split(b[i],c,"+");
		donor=c[1];
		acceptor=c[2];
		revDonor=revComp(donor);
		revAcceptor=revComp(acceptor);	

		regexDonor="^"donor;
		regexAcceptor=acceptor"$";
		regexRevDonor=revDonor"$";
		regexRevAcceptor="^"revAcceptor;

		if ((seq1 ~ regexDonor)&&(seq2 ~ regexAcceptor))
		{
			juncDonor=donor;
			juncAcceptor=acceptor;
			juncStr1=mapStr1;
			juncStr2=mapStr2;		
			rev=0; # Blocks in biological order: first block donor and second one acceptor. 
			break;
		}
		else if ((seq2 ~ regexRevDonor)&&(seq1 ~ regexRevAcceptor))
		{
			juncDonor=donor;
			juncAcceptor=acceptor;
			juncStr1=(mapStr1=="+" ? "-" : "+");
			juncStr2=(mapStr2=="+" ? "-" : "+");
			rev=1; # Blocks in not biological order: first block acceptor and second one donor. 
			break;
		}
	}
	
	if ((juncStr1=="")||(juncStr2=="")) # If no consensus splice sites (css). However this case should not exist, since the rna-mapper requires css to split-map the reads
	{
		print "[ERROR] Consensus splice sites were not found." > "/dev/stderr"
        exit 1
	}
}

function makeChimJunc(chr1, breakpoint1, juncStr1, chr2, breakpoint2, juncStr2, rev)
{
	if (rev=="0")
	{	
		chimJunc=chr1"_"breakpoint1"_"juncStr1":"chr2"_"breakpoint2"_"juncStr2;	
	}
	else # Reverse to write first the donnor block and then the acceptor
	{	
		chimJunc=chr2"_"breakpoint2"_"juncStr2":"chr1"_"breakpoint1"_"juncStr1;		
	}
	return chimJunc
}

function findMaxBegEnd(beg, end, chimJunc)
{
	split(chimJunc,a,":");
	split(a[1],b1,"_");
	split(a[2],b2,"_");
	str1=b1[3];
	str2=b2[3];
	
	if ((maxBeg[chimJunc]=="")||(maxEnd[chimJunc]==""))
	{
		maxBeg[chimJunc]=beg;
		maxEnd[chimJunc]=end;
	}
	else
	{
		if (str1=="+")
		{
			if (str2=="+") # junc +/+
			{
				maxBeg[chimJunc]=(maxBeg[chimJunc] > beg ? beg : maxBeg[chimJunc]);
				maxEnd[chimJunc]=(maxEnd[chimJunc] < end ? end : maxEnd[chimJunc]);
			}
			else # junc +/-
			{
				maxBeg[chimJunc]=(maxBeg[chimJunc] > beg ? beg : maxBeg[chimJunc]);
				maxEnd[chimJunc]=(maxEnd[chimJunc] > end ? end : maxEnd[chimJunc]);
			}
		}		
		else
		{
			if (str2=="+") # junc -/+
			{
				maxBeg[chimJunc]=(maxBeg[chimJunc] < beg ? beg : maxBeg[chimJunc]);
				maxEnd[chimJunc]=(maxEnd[chimJunc] < end ? end : maxEnd[chimJunc]);
			}
			else # junc -/-
			{
				maxBeg[chimJunc]=(maxBeg[chimJunc] < beg ? beg : maxBeg[chimJunc]);
				maxEnd[chimJunc]=(maxEnd[chimJunc] > end ? end : maxEnd[chimJunc]);
			}
		}
	}
}

{ 
	blocks=$1;
   	split(blocks,b,":"); 
  	split(b[1],b1,"_"); 
	split(b[2],b2,"_"); 	
	chr1=b1[1];
	chr2=b2[1];
	mapStr1=b1[4];
	mapStr2=b2[4];

	if (mapStr1=="+")
	{
		if (mapStr2=="+") # 2 Blocks mapped in +/+
		{
   			# Chimeric junction breakpoints
   			breakpoint1=b1[3];
			breakpoint2=b2[2];
			
			# 8-mers containing the donor and acceptor splice sites 
			seq1=$4;
			seq2=$5;
		
			# Find the splice sites in the 8-mers and determine the strand of the blocks based on the splice sites sequences
			findStrand(seq1,seq2,mapStr1,mapStr2,CSS);
			
			# Make the chimeric junction
			makeChimJunc(chr1, breakpoint1, juncStr1, chr2, breakpoint2, juncStr2, rev);
			strandedness
			# Determine the maximum external coordinates of the junction 		
			beg=(rev=="0" ? b1[2] : b2[3]);
			end=(rev=="0" ? b2[3] : b1[2]);
			
			findMaxBegEnd(beg, end, chimJunc);
		}
		else # 2 Blocks mapped in +/-
		{
			# Chimeric junction breakpoints
			breakpoint1=b1[3];
			breakpoint2=b2[3];
		
			# 8-mers containing the donor and acceptor splice sites 
			seq1=$4;
			seq2=$6;
			
			# Find the splice sites in the 8-mers and determine the strand of the blocks based on the splice sites sequences
			findStrand(seq1,seq2,mapStr1,mapStr2,CSS);
			
			# Make the chimeric junction
			makeChimJunc(chr1, breakpoint1, juncStr1, chr2, breakpoint2, juncStr2, rev);
			
			# Determine the maximum external coordinates of the junction 		
			beg=(rev=="0" ? b1[2] : b2[2]);
			end=(rev=="0" ? b2[2] : b1[2]);
			
			findMaxBegEnd(beg, end, chimJunc);
		}
	}
	else
	{
		if (mapStr2=="+") # 2 Blocks mapped in -/+
		{
			# Chimeric junction breakpoints
			breakpoint1=b1[2];
			breakpoint2=b2[2];
		
			# 8-mers containing the donor and acceptor splice sites 
			seq1=$3;
			seq2=$5;
					
			# Find the splice sites in the 8-mers and determine the strand of the blocks based on the splice sites sequences
			findStrand(seq1,seq2,mapStr1,mapStr2,CSS);
			
			# Make the chimeric junction
			makeChimJunc(chr1, breakpoint1, juncStr1, chr2, breakpoint2, juncStr2, rev);
			
			# Determine the maximum external coordinates of the junction 		
			beg=(rev=="0" ? b1[3] : b2[3]);
			end=(rev=="0" ? b2[3] : b1[3]);
			
			findMaxBegEnd(beg, end, chimJunc);			
		}
		else # 2 Blocks mapped in -/-
		{
			# Chimeric junction breakpoints
			breakpoint1=b1[2];
			breakpoint2=b2[3];
		
			# 8-mers containing the donor and acceptor splice sites 
			seq1=$3;
			seq2=$6;
					
			# Find the splice sites in the 8-mers and determine the strand of the blocks based on the splice sites sequences
			findStrand(seq1,seq2,mapStr1,mapStr2,CSS);
			
			# Make the chimeric junction
			makeChimJunc(chr1, breakpoint1, juncStr1, chr2, breakpoint2, juncStr2, rev);
			
			# Determine the maximum external coordinates of the junction 		
			beg=(rev=="0" ? b1[3] : b2[2]);
			end=(rev=="0" ? b2[2] : b1[3]);
			
			findMaxBegEnd(beg, end, chimJunc);
		}
	}	
	if ((strandedness!="1")||(rev!="1")) # Filter out junction if it is stranded data and the read split-map with the complementary
									 # reverse of the consensus splice sites. This case is a gem mapping artefact.
	{												
		# Save the donor and the acceptor associated to the chimeric junction into a dictionary
		juncDonors[chimJunc]=juncDonor;
		juncAcceptors[chimJunc]=juncAcceptor;
	
		# Count the number of staggered reads and the total number of reads supporting the junction
		nbStag[chimJunc]++;
		nbTotal[chimJunc]=nbTotal[chimJunc]+$2;
	}
}		
	
END{
	for (chimJunc in nbStag)
	{
		print chimJunc, nbStag[chimJunc], nbTotal[chimJunc], maxBeg[chimJunc], maxEnd[chimJunc], juncDonors[chimJunc], juncAcceptors[chimJunc];
	}
}


