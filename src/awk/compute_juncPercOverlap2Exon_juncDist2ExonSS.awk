#!/usr/bin/env awk

# *****************************************************************************
	
#	compute_JuncPercOverlap_JuncDistance2ExonSS.awk

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

# Takes as input a file resulting of intersecting splice junction side sequences with annotated exons through bedtools intersectBed. Computes for those splice junction sides overlapping exons, the percentage overlapping the exon and the distance between the splice junction donor/acceptor to the exon boundary.

## input
# chr3	ChimPipe	spliceJunc1	182637098	182637125	.	+	.	JuncId: "chr3_182637125_+:chr3_182637126_+";	chr3	HAVANA	exon	182635811	182639423	.	+	.	gene_name "ATP11B"; gene_id "ENSG00000058063.11"; gene_type "protein_coding"; transcript_id "ENST00000323116.5";28
# chr3	ChimPipe	spliceJunc2	182637126	182637150	.	+	.	JuncId: "chr3_182637125_+:chr3_182637126_+";	chr3	HAVANA	exon	182635811	182639423	.	+	.	gene_name "ATP11B"; gene_id "ENSG00000058063.11"; gene_type "protein_coding"; transcript_id "ENST00000323116.5";25
# chr4	ChimPipe	spliceJunc1	71629354	71629386	.	+	.	JuncId: "chr4_71629386_+:chr4_71630192_+";	chr4	HAVANA	exon	71629269	71629386	.	+	.	gene_name "RUFY3"; gene_id "ENSG00000018189.8"; gene_type "protein_coding"; transcript_id "ENST00000503876.1";	33
# chr4	ChimPipe	spliceJunc2	71630192	71630235	.	+	.	JuncId: "chr4_71629386_+:chr4_71630192_+";	chr4	HAVANA	exon	71630192	71630293	.	+	.	gene_name "RUFY3"; gene_id "ENSG00000018189.8"; gene_type "protein_coding"; transcript_id "ENST00000503876.1";	44
# chr20	ChimPipe	spliceJunc1	6855536	6855539	.	+	.	JuncId: "chr20_6855539_+:chr20_24943838_+";	.	.	.	-1	-1	-1	.	.	.	0

## output
# chr3	ChimPipe	spliceJunc1	182637098	182637125	.	+	.	JuncId: "chr3_182637125_+:chr3_182637126_+";	chr3	HAVANA	exon	182635811	182639423	.	+	.	gene_name "ATP11B"; gene_id "ENSG00000058063.11"; gene_type "protein_coding"; transcript_id "ENST00000323116.5";28 100 2298
# chr3	ChimPipe	spliceJunc2	182637126	182637150	.	+	.	JuncId: "chr3_182637125_+:chr3_182637126_+";	chr3	HAVANA	exon	182635811	182639423	.	+	.	gene_name "ATP11B"; gene_id "ENSG00000058063.11"; gene_type "protein_coding"; transcript_id "ENST00000323116.5";25 100 1315
# chr4	ChimPipe	spliceJunc1	71629354	71629386	.	+	.	JuncId: "chr4_71629386_+:chr4_71630192_+";	chr4	HAVANA	exon	71629269	71629386	.	+	.	gene_name "RUFY3"; gene_id "ENSG00000018189.8"; gene_type "protein_coding"; transcript_id "ENST00000503876.1";	33 100 0
# chr4	ChimPipe	spliceJunc2	71630192	71630235	.	+	.	JuncId: "chr4_71629386_+:chr4_71630192_+";	chr4	HAVANA	exon	71630192	71630293	.	+	.	gene_name "RUFY3"; gene_id "ENSG00000018189.8"; gene_type "protein_coding"; transcript_id "ENST00000503876.1";	44 100 0
# chr20	ChimPipe	spliceJunc1	6855536	6855539	.	+	.	JuncId: "chr20_6855539_+:chr20_24943838_+";	.	.	.	-1	-1	-1	.	.	.	0 -1 -1


# Usage example:
################
# awk -f compute_JuncPercOverlap_JuncDistance2ExonSS.awk spliceJunctions_2parts_beg_end_intersected.txt

## A) Splice junction side does not overlap any annotated exon
$NF==0{
	print $0, "-1", "-1";
}

## B) Splice junction side overlap annotated exon
$NF!="0"{
	
	# B.1 Compute percentage of the splice junction side sequence that overlap the annotated exon
	blockLength=($5-$4+1); 
	nbOverlappingBases=$NF; 
	percOverlap=(nbOverlappingBases/blockLength)*100; 
	
	# B.2 Compute distance between the splice junction donor/acceptor to the exon boundary
	juncStrand=$7; 

	if ((($3=="spliceJunc1")&&(juncStrand=="+"))||(($3=="spliceJunc2")&&(juncStrand=="-")))
	{
		juncCoord=$5; 
		exonBoundaryCoord=$15;
	}
	else
	{
		juncCoord=$4; 
		exonBoundaryCoord=$14;
	}; 

	dist=juncCoord-exonBoundaryCoord;  
	exonBoundDist=(dist < 0 ? -dist : dist); 

	# B.3 Print output
	print $0, percOverlap, exonBoundDist;
}



