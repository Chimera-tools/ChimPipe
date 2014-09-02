# 31/07/2014
#

# awk -v CSS='GT+AG,GC+AG,ATATC+A.,GTATC+AT' -f make_chimjunctions.awk input

# input
# chr10_100018895_100018913_-:chr3_58279376_58279407_+	1	GGAGCCAC	ATGTGGGG	AGGTGGAG	GGAAAGAC		
# chr10_100143466_100143480_+:chr3_183901336_183901371_+	1	GTACTAAC	CTGGGAAG	TGAACTAC	GGCCGCAG		
# chr10_100154777_100154810_+:chr19_36214052_36214068_+	1	CGCAGTTC	GTGATGTG	ACGGCAAG	ATGTGGAC		
# chr10_100171018_100171033_+:chr7_72419139_72419173_-	2	GCTTCCAT	GTGATGGT	GCTGAGCT	CTGCGCCT		


# Output
# chr11_60781486_+:chr11_120100337_+ 1 1 60781475 120100375 GT AG
# chr16_30537851_-:chr15_40629945_+ 1 1 30537885 40629960 AC AG
# chr13_64414915_-:chr18_74690910_- 1 3 64414930 74690876 AC CT
# chr11_71714505_+:chr22_50928336_- 1 1 71714472 50928320 GT CT


# Position of the splice sites regarding on the strand the blocks map. 
######################################################################
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

#		print "----- forward -----";
#		print seq1, regexDonor;
#		print seq2, regexAcceptor;
#		print "----- rev -----";
#		print seq2, regexRevDonor;
#		print seq1, regexRevAcceptor;

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
	
	if ((juncStr1=="")||(juncStr2==""))
	{
		juncStr1=".";
		juncStr2=".";
		juncDonor="NA";
		juncAcceptor="NA";
	}
}

function makeChimJunc(chr1, breakpoint1, juncStr1, chr2, breakpoint2, juncStr2, rev)
{
	if (rev=="0")
	{	
		chimJunc=chr1"_"breakpoint1"_"juncStr1":"chr2"_"breakpoint2"_"juncStr2;	
	}
	else
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
#	print mapStr1, mapStr2;
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
	if ((stranded!="1")||(rev!="1"))
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


