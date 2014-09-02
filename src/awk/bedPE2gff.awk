
## Input: BEDPED
#chr1 752991 753018 chr1 754245 754250 HWI-ST985:73:C08BWACXX:8:2208:2017:40383/1  1000 + +

## Output: GFF
# chr1 hts alblock 752992 753018 . + . name: HWI-ST985:73:C08BWACXX:8:2208:2017:40383/1
# chr1 hts alblock 754246 754250 . + . name: HWI-ST985:73:C08BWACXX:8:2208:2017:40383/1

BEGIN{OFS="\t"}
{
	block1=$1"\thts\talblock\t"($2+1)"\t"$3"\t.\t"$9"\t.\tname: \""$7"\"\;";
	block2=$1"\thts\talblock\t"($5+1)"\t"$6"\t.\t"$10"\t.\tname: \""$7"\"\;";
	
	print block1;
	print block2;
}
