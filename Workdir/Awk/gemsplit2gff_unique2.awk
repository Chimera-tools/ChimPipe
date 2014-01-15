# ~/Awk/gemsplit2gff_unique2.awk

# very different from ~/Awk/gemsplit2gff_unique.awk which was done by Vincent Lacroix on the first version of the gem in 2008
# split mapper. Here completely rewritten to deal with the new convention on writting the split mappings

# takes as input a split mapping gem file from gem-rna-mapper (June 2013) and outputs the unique 2 block split mappings
# in gff format = each block on a different gff row (consecutive for a given split mapping) and with the read id in column 10

# example of input
##################
# SINATRA_0006:1:1:7430:930#0/1	NCCTTCTCTTCGCTCCTGGTGTAAGGTATGGTACATAAGAGTCCAATGCTATTTGCGCAAGTGCTAGGGTAACGAG	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	0:0:0:0:1	chr1:+:80324753:15::chr12:+:112677685:61
# SINATRA_0006:1:1:3790:968#0/1	NAACTCATCATAGTGTTCCTGCATCTCCACATCGCTCACGGCACAGTGTGAGCCGTCAGCCGTCTGTGCACTGTTT	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	0:0:1	chr21:+:44515810:44::chr21:+:44521476:32

# resulting output
##################
# chr1	hts	alblock	80324738	80324753	.	+	.	name: "SINATRA_0006:1:1:7430:930#0/1";
# chr12	hts	alblock	112677685	112677746	.	+	.	name: "SINATRA_0006:1:1:7430:930#0/1";
# chr21	hts	alblock	44515767	44515810	.	+	.	name: "SINATRA_0006:1:1:3790:968#0/1";
# chr21	hts	alblock	44521476	44521507	.	+	.	name: "SINATRA_0006:1:1:3790:968#0/1";


$NF!="-"{
    n=split($NF,a,","); 
    if(n==1)
    {
	n2=split($NF,b,"::"); 
	if(n2==2)
	{
	    split(b[1],b1,":"); 
	    split(b[2],b2,":"); 
	    print b1[1], "hts", "alblock", b1[3]-b1[4]+1, b1[3],  ".", b1[2], ".", "name:", "\""$1"\"\;"; 
	    print b2[1], "hts", "alblock", b2[3], b2[3]+b2[4]-1,  ".", b2[2], ".", "name:", "\""$1"\"\;";
	}
    }
}