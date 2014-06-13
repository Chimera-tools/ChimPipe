# ~/Awk/bed12fields2gff.awk
# usually used to process alignments of rna sequencing data
# improved on 07/02/2013 to be able to write the blocks in reverse order 
# this is to deal with bed12 obtained from bam where the convention is to write the split mappings 
# of the - strand in the order relative to the + strand (unlike in gem split mapper format)

# awk 'NR>=2' GM12878rnaseqr1.all.bed | awk -v rev=0 -f ~/Awk/bed12fields2gff.awk | gff2gff

# awk 'NR>=2' GM12878rnaseqr1.all.bed 
# chr1    752991  754250  7488914 1000    +       0       0       0,0,255 2       27,5    0,1254
 
# output
# chr1    hts     alblock 752992  753018  .       +       .       name: 7488914
# chr1    hts     alblock 754246  754250  .       +       .       name: 7488914

{
    split($11,sz,","); 
    split($12,st,",");
    for(i=1; i<=$10; i++)
    {
	line[i]=$1"\thts\talblock\t"($2+st[i]+1)"\t"($2+st[i]+sz[i])"\t.\t"$6"\t.\tname: \""$4"\"\;";
    }
    if ((rev=="")||(rev==0)||($6=="+")) 
    {
	for(i=1; i<=$10; i++)
	{
	    print line[i];
	}
    }
    else
    {
	for(i=$10; i>=1; i--)
	{
	    print line[i];
	}
    }
}
