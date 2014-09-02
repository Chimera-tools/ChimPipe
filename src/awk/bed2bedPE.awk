#!/usr/bin/env awk

# Convert into BED paired-end format a bed12 file containing split-mapped reads in two blocks. If input variable rev set to different than 0, 
# reverse the blocks when reads split-mapping in -/-

## Input: BED12
#chr22 1000 5000 HWI-ST985:73:C08BWACXX:8:2208:2017:40383/1 960 + 1000 5000 0 2 567,488, 0,3512

## Output: BEDPED
#chr22 1000 1567 chr22 4512 5000 HWI-ST985:73:C08BWACXX:8:2208:2017:40383/1 960 + +

# Usage example:
# awk -v rev="1" -f bed2bedPE.awk

BEGIN{OFS="\t"}
{
	split($11,sz,",");
	split($12,st,",");
	begA=$2+st[1];
	endA=begA+sz[1];
	begB=$2+st[2];
	endB=begB+sz[2];
	if ((rev=="")||(rev==0)||($6=="+")) 
    {
		print $1, begA, endA, $1, begB, endB, $4, $5, $6, $6;
	}
	else
	{
		print $1, begB, endB, $1, begA, endA, $4, $5, $6, $6;
	}
}



