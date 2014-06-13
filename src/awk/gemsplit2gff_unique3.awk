# ~/Awk/gemsplit2gff_unique3.awk

# takes as input a split mapping gem file from gem-rna-mapper (June 2013) and outputs the unique 2 block split mappings
# in gff format = each block on a different gff row (consecutive for a given split mapping) and with the read id in column 10
# but different from ~/Awk/gemsplit2gff_unique2.awk to be correct and use Paolo's code
# cat MY_FILE.gem | awk -F '\t' '{split($5,s,","); len=length(s); for (i=1;i<=len;++i) {split(s[i],ss,"::"); split(ss[1],ss1,":"); split(ss[2],ss2,":"); print ss1[1]"\thts\talblock\t"ss1[3]"\t"(ss1[3]+ss1[4]-1)"\t.\t"ss1[2]"\t.\tname:\""$1"\";"; print ss2[1]"\thts\talblock\t"ss2[3]"\t"(ss2[3]+ss2[4]-1)"\t.\t"ss2[2]"\t.\tname:\""$1"\";"}}'

# example of input
##################
# SINATRA_0006:1:1:7430:930#0/1	NCCTTCTCTTCGCTCCTGGTGTAAGGTATGGTACATAAGAGTCCAATGCTATTTGCGCAAGTGCTAGGGTAACGAG	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	0:0:0:0:1	chr1:+:80324753:15::chr12:+:112677685:61
# SINATRA_0006:1:1:3790:968#0/1	NAACTCATCATAGTGTTCCTGCATCTCCACATCGCTCACGGCACAGTGTGAGCCGTCAGCCGTCTGTGCACTGTTT	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	0:0:1	chr21:+:44515810:44::chr21:+:44521476:32



$NF!="-"{
    n=split($NF,a,","); 
    if(n==1)
    {
	n2=split($NF,b,"::"); 
	if(n2==2)
	{
	    split(b[1],b1,":"); 
	    split(b[2],b2,":"); 
	    print b1[1], "hts", "alblock", b1[3], b1[3]+b1[4]-1,  ".", b1[2], ".", "name:", "\""$1"\"\;"; 
	    print b2[1], "hts", "alblock", b2[3], b2[3]+b2[4]-1,  ".", b2[2], ".", "name:", "\""$1"\"\;";
	}
    }
}
