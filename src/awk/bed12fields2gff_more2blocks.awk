#!/bin/bash
# Adapted from: ~sdjebali/Awk/bed12fields2gff.awk to deal with more than 2 blocks.

{
    split($11,sz,","); 
    split($12,st,",");
    for(i=1; i<=$10; i++)
    {
	line[i]=$1"\thts\talblock\t"($2+st[i]+1)"\t"($2+st[i]+sz[i])"\t.\t"$6"\t.\tname: \""$4"\"\;";
    }
    if(rev==""||rev==0)
    {
	for(i=1; i<$10; i++)
	{
	    print line[i];
	    print line[i+1];
	}
    }
    else
    {
	for(i=($10-1); i>=1; i--)
	{
	    print line[i+1];
	    print line[i];
	}
    }
}
