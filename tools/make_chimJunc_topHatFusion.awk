#!/usr/bin/env awk

# Description
###############
# Takes as input TopHat output file for chimeric junction (--fusion-search option has to be enabled) and produces chimeric junctions in ChimPipe's format plus their support

### Input: TopHat-fusion output file for chimeric junctions (fusions.out file)
# chr20-chr17 49411707 59445685 ff 106 116 167 0 37 36 0.569598 11 25 38 49 63 CAGCGGGG........ GGTTAGGT........ 106 106 106 106 106 106 106 106..... 106 106 106 106 106 106 106 106....... -6:1 11:0 16:-3 18:1 14:6 14:6 15:7 23:0 5:21 31:5 18:19 36:-1 ...

# Fusion between 49411707th base on chromosome 20 and 59445685th base on chromosome 17.
# ff is the orientations of the two chromosomes - both chromosomes are in forwarding direction, like (chr 20) -----> -----> (chr 17).
# 106 is the number of reads that span the fusion (nbSupporingReads), 116 is the number of mate pairs that support the fusion (nbSupportingPairs), 167 is the number of mate pairs that support the fusion and whose one end spans the fusion.
# 0 is the number of reads that contradict the fusion by mapping to only one of the chromosomes 20 and 17.
# 37 and 36 are the number of bases on the left and right sides of a fusion, respectively, covered by spanning reads.
# The nucleotide sequences correspond to the two 50-bp contigs around a fusion point on chromosome 20. 
# The next fields group is the depth coverage by spanning reads on chromosome 20 from 49411707th base to 49411658th base. The sixth row is depth coverage by spanning reads on chromosome 17 from 59445685th base to 59445734th base.
# The last fields group is the distances (distance1:distance2) between a mate pair and the fusion, i.e., distance1 bewtween the left end of a pair and the left side of the fusion and distance2 between the right end of a pair and the right side of the fusion. 

### Output: Chimeric junctions in ChimPipe's format + number of reads supporting each junction
# chimJunc nbSupportingReads	nbSupportingPairs	nbSupportingReadsConsistentPair	nbContradictoryReads
# chr4_76807216_+:chr2_191402746_- 10	5	8	0	

# Usage:
########
# awk -f make_chimJunc_topHatFusion_output.awk topHatFusion_output_chimeras


BEGIN{
    OFS="\t";
    print "chimJunc", "nbSupportingReads", "nbSupportingPairs", "nbSupportingReadsConsistentPair", "nbContradictoryReads";
}

! /^#/{
	split($1,chr,"-"); 
	split($4,strand,""); 
	
	### 5 prime junction boundary
	chr5prime=chr[1];
	coord5prime=$2;
		
	# Convert topHat-fusion strand information (f=="+"; r=="-") into +,-
	strand5prime=(strand[1]=="f" ? "+" : "-");

	### 3 prime junction boundary
	chr3prime=chr[2];
	coord3prime=$3;
	
	# Convert topHat-fusion strand information (f=="+"; r=="-") into +,-
	strand3prime=(strand[2]=="f" ? "+" : "-");
	
	### Make chimJunc
	chimJunc=chr5prime"_"coord5prime"_"strand5prime":"chr3prime"_"coord3prime"_"strand3prime;
	
	### Print output
	
	print chimJunc, $5, $6, $7, $8; 
}

