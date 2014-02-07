# ~/Awk/gemsplit2gff_unique4.awk

# takes as input a split mapping gem file from gem-rna-mapper (June 2013) and outputs the unique 2 block split mappings
# in gff format = each block on a different gff row (consecutive for a given split mapping) and with the read id in column 10
# but different from ~brodriguez/Chimeras_project/Chimeras_detection_pipeline/Chimera_mapping/Workdir/Awk/gemsplit2gff_unique3.awk since it is able to write the blocks in reverse order when the split-mappings are in the same chromosome, same strand and in the minus strand. This is to generate a gff file with the same convention than the
# one generated from a bam file, since bam and gem has different conventions, while the bam provides the alignment in genomic order gem does it in the biological one (from 5' to 3').

#  awk -v rev=0 -f ~/Awk/gemsplit2gff_unique4.awk

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
	    	if ((b1[1]!=b2[1]) || (b1[2]!=b2[2]) || (b1[2]!="-") || (rev=="") || (rev==0))
	    	{
	    		print b1[1], "hts", "alblock", b1[3], b1[3]+b1[4]-1,  ".", b1[2], ".", "name:", "\""$1"\"\;"; 
	   	 		print b2[1], "hts", "alblock", b2[3], b2[3]+b2[4]-1,  ".", b2[2], ".", "name:", "\""$1"\"\;";
			}
			else 
			{ 
	   	 		print b2[1], "hts", "alblock", b2[3], b2[3]+b2[4]-1,  ".", b2[2], ".", "name:", "\""$1"\"\;";
				print b1[1], "hts", "alblock", b1[3], b1[3]+b1[4]-1,  ".", b1[2], ".", "name:", "\""$1"\"\;";
			}	
		}
	}
}
