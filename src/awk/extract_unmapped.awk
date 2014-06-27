#!/usr/bin/env awk

# Take as input a map file containing mapped reads with the gem mapper, extract the reads that do not map with a number of mismatches lower than 6 and produce a map file with them

# Input (mapped reads in map format)
####################################
#SRR201779.2_PATHBIO-SOLEXA2_30TUEAAXX:3:1:966:1023_length=53#0/2	TGATCCTTAGGAGTGTTTCCTCGCTTAAGGGAAATTAGGTCCCATGATGATGA	BBBBCCCBBBBBCCCBCCCCCCCBBBBBBBBBBBBBBBBBBBBBBBBBBBB?<	0	-
#SRR201779.3_PATHBIO-SOLEXA2_30TUEAAXX:3:1:929:989_length=53#0	TGATACAATGTATGAAATATTTGTGTCTGTACAGGGTTAAAGTGGCTGAAGTT TACTACATCCGAGCAGTGTCTGATAGATGGTTGGGTGCTGAGGCAGTATGTAT	BBBBBBBB8BBBBBBBBB?BBBBBB<BBBBBBBBA@@>@BB@@AB<9969<7? BCBBBBBBBBBBBBBBBBBBBBBBBBBA@@B@A@BBBBB???B@=???<??<?	1	chr6:+:101086580:53::chr6:-:101090507:53:::32736

# Output (unmapped reads)
####################################
#SRR201779.32_PATHBIO-SOLEXA2_30TUEAAXX:3:1:1097:958_length=53#0/1	AAAAAGGACAAGACTGTATGACTTCTATGTTTGCCCCGGTCATACTGTACCAA	BBB6(9BBBBBCBBBBBBBBBCCBCBBBBBBBBBBBBBBBBBBBB<9<BBBB@	0	-
#SRR201779.38_PATHBIO-SOLEXA2_30TUEAAXX:3:1:924:811_length=53#0	AGAGGTTAGGATAATGGAAGTAGAAAGTCCTGTTTGGAAAAGGGAGATCGGAA CCTTTTCCAAACAGGACTTTCTACTTCCATTATCCTAACCTCTCAGATCGGAA	BBBBBBBBBBBBBBBBBBB?<<@@BBAAABA>@AA>ABBB@?;@@A>@<<?;< CCCCCCCBBBCBBBBBBBBBBBBBBBBBBBBBBBBB@;@BBBBA@@@ABBA?9	0:0:0:0:0:0:0:0:0:1	chr1:+:10106499:45>1-1>2-1>1-2::chr1:-:10106488:43>1+2CCA1T1T1:::32640

{
	if (NF==5)
	{
		m=$4;
	}
	else
	{
		m=$6
	}; 
	gsub(/+/, ":", m);  # In some cases the gem alignment format has a + in place of a colon in the match summary. I replace it by a colon to 
	split(m,a,":");     # only have colons and avoid problems
	n=1; 
	discard=0; 
	while (n < 7 )
	{
		if (a[n]!=0)
		{
			discard=1;
		}
		n++
	}
	if (discard!=1)
	{
		print $0
	}
}


