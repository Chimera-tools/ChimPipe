#!/usr/bin/env awk

# *****************************************************************************
	
#	select_overlappingGene_contiguousMapping.awk
	
#	This file is part of the ChimPipe pipeline 

#	Copyright (c) 2014 Bernardo Rodríguez-Martín 
#					   Emilio Palumbo 
#					   Sarah djebali 
	
#	Computational Biology of RNA Processing group
#	Department of Bioinformatics and Genomics
#	Centre for Genomic Regulation (CRG)
					   
#	Github repository - https://github.com/Chimera-tools/ChimPipe
	
#	Documentation - https://chimpipe.readthedocs.org/

#	Contact - chimpipe.pipeline@gmail.com
	
#	Licenced under the GNU General Public License 3.0 license.
#******************************************************************************

# Description
###############

# 

## input
# chr1    ChimPipe        alBlock 10525   10577   1       -       .       ReadName: "SRR201779.5012779_PATHBIO-SOLEXA2_30TUEAAXX:3:78:1387:1378_length=53#0/2";   .       .       .       -1  -1       .       .       .       .       0
# chr1    ChimPipe        alBlock 11608   11660   1       +       .       ReadName: "SRR201779.1488862_PATHBIO-SOLEXA2_30TUEAAXX:3:22:1724:187_length=53#0/2";    .       .       .       -1  -1       .       .       .       .       0
# chr1    ChimPipe        alBlock 12118   12170   1       -       .       ReadName: "SRR201779.1488862_PATHBIO-SOLEXA2_30TUEAAXX:3:22:1724:187_length=53#0/1";    chr1    HAVANA  exon    1186912227   .       +       .       gene_id "ENSG00000223972.4"; gene_name "DDX11L1"; gene_type "pseudogene"; transcript_id "ENST00000456328.2";    53
# chr1    ChimPipe        alBlock 12118   12170   1       -       .       ReadName: "SRR201779.1488862_PATHBIO-SOLEXA2_30TUEAAXX:3:22:1724:187_length=53#0/1";    chr1    ENSEMBL exon    1187212227   .       +       .       gene_id "ENSG00000223972.4"; gene_name "DDX11L1"; gene_type "pseudogene"; transcript_id "ENST00000515242.2";    53

## output
# SRR201779.5012779_PATHBIO-SOLEXA2_30TUEAAXX:3:78:1387:1378_length=53#0/2        chr1_10525_10577_-      unannotated
# SRR201779.1488862_PATHBIO-SOLEXA2_30TUEAAXX:3:22:1724:187_length=53#0/2 chr1_11608_11660_+      unannotated
# SRR201779.1488862_PATHBIO-SOLEXA2_30TUEAAXX:3:22:1724:187_length=53#0/1 chr1_12118_12170_-      ENSG00000223972.4
# SRR201779.3108351_PATHBIO-SOLEXA2_30TUEAAXX:3:48:779:1690_length=53#0/2 chr1_17291_17343_+      ENSG00000227232.4


# Usage example:
################
# awk -f 


{	 	
	# Set the mapping block id as the read name
	split($10,a,"\""); 
	readId=a[2];
	
	# Read mapping coordinates
	mapCoord=$1"_"$4"_"$5"_"$7;
	
	
	# A) Read do not overlapping any annotated exon
	if ($NF == "0")
	{
		gnId="unannotated"; 
	}
		
	# B) Read overlapping an annotated exon
	else
	{
		# Parse exon atributes list and select its gene id
		for(i=19; i<=(NF); i++)
		{
			if ($i=="gene_id")
			{
				split($(i+1),id,"\""); 
				gnId=id[2]; 
			}
		}		
	};

	# Print output:
	row=readId"\t"mapCoord"\t"gnId;
	
	print row;
}
