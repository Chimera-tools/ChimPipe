#!/usr/bin/env awk

# Description
###############
# Takes as input FusionMap output file for chimeric junction (split-read analysis mode) and produces chimeric junctions in ChimPipe's format plus their support

### Input: FusionMap output file for chimeric junctions (*.FusionReport.txt file)
# FusionID	DatasetP2_SimulatedReads.UniqueCuttingPositionCount	DatasetP2_SimulatedReads.SeedCount	DatasetP2_SimulatedReads.RescuedCount	Strand	Chromosome1	Position1	Chromosome2	Position2	KnownGene1	KnownTranscript1	KnownExonNumber1	KnownTranscriptStrand1	KnownGene2	KnownTranscript2	KnownExonNumber2	KnownTranscriptStrand2	FusionJunctionSequence	FusionGene	SplicePattern	SplicePatternClass	FrameShift	FrameShiftClass	Distance	OnExonBoundary	Filter
# FUS_11184689_446432191(++)	27	32	6	--	11	61646784	1	11184690	FADS3	NM_021727_chr11	3	-	MTOR	NM_004958_chr1	47	-	AGACAAATGTAGGAAAAAACCAGAAGACTTCTCAAATTGTTGCCATTTCAGGGTTTCTGAATACCTGAGGTTTTTCCGAAGAGATGTTGGGTCATTGGCCAGAAGGGTGTTAACCAGGCCGAAGAGCTGCATCACACGCTCATCCTGGCGCAGATCTTCATGGCCTTTTAGAAGGAAAACAAACTCATGTCCGTTGCTGC@ctgagagatggccaggatgaaggcggccagggcactgggcacccagccaggacccaggaggtagataaggagccaggccagcacctccatggccaggatgtggcccagtaggaaagcaaagaaggtgggactggcatcaaacagcttcatgtcctcggctgcctggtgcagggctcggaagtcctcgaccagctgcgcct	FADS3->MTOR	GT-AG	CanonicalPattern[Major]	2->1	FrameShift	-1	Both

#### Meaning different type of supporting counts:
# source: http://www.arrayserver.com/wiki/index.php?title=FusionMap

## A) Unique cutting position count: Indicates the number of unique cutting positions of a specified FusionID. It is expected that a real fusion would have more than 1 unique cutting positions.

# For instance, if a fusion occurs across two genes as follows:

# Read1 AAATTTCCCCCGTTTTTGACCCCCCCCC

# Read1: Gene1 AAATTTCCCCCG Break Gene2: TTTTTGACCCCCCCCC

# It would be expected that the start position of the read in Gene1 would differ if the fusion is real.

# Read2 AAAAAATTTCCCCCGTTTTTGACCCCCC

# Read2: Gene1 AAAAAATTTCCCCCG "Break" Gene2 TTTTTGACCCCCC

# In the above example, there are two different unique cutting positions.

# Unique cutting position is the same as unique mapping location of fusion reads, and usually indicate unique read sequences

## B) SeedCount: number of seed supporting reads. A seed read is a read that has at least the MinimalFusionAlignmentLength for each gene, 

## C) RescuedCount: number of rescued supporting reads. A rescued read is a read that does not have at least the MinimalFusionAlignmentLength for each gene

# Usage:
########
# awk -f make_chimJunc_fusionMap_output.awk fusionMap_output_chimeras


BEGIN{
    OFS="\t";
    print "chimJunc", "uniqueCuttingPositionCount", "seedCount", "rescuedCount";
}

NR>1{
	split($5,strands,"");
	strandA=strands[1];
	strandB=strands[2];
	chimJunc="chr"$6"_"$7"_"strandA":chr"$8"_"$9"_"strandB;
	uniqueCuttingPos=$2;
	seedCount=$3;
	rescuedCount=$4;
	
	print chimJunc, uniqueCuttingPos, seedCount, rescuedCount;
}

