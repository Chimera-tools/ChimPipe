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

# Takes as input the resulting output file of the intersection between the alignments and the annotated exons. Colapse the alignment coordinates in a single column, select the gene id and flag as unannotated those aligment do not overlapping any exon. 

## input
# chr1    14551   14601   SRR064286.5709675_HWI-EAS418:3:7:315:351_length=50#0/1  1	+	chr1    14362   14829   exon:ENSG00000227232.4;WASH7P;pseudogene        .	-	50
# chr1    14681   14731   SRR064286.5709675_HWI-EAS418:3:7:315:351_length=50#0/2  1	-	chr1    14362   14829   exon:ENSG00000227232.4;WASH7P;pseudogene        .	-	50
# chr1    14912   14962   SRR064286.7043546_HWI-EAS418:3:53:1606:39_length=50#0/1 1	+	.	-1	-1	.	-1	.	0
# chr1    14928   14978   SRR064286.7043546_HWI-EAS418:3:53:1606:39_length=50#0/2 1	-	chr1    14969   15038   exon:ENSG00000227232.4;WASH7P;pseudogene        .	-	9

## output
# SRR064286.5709675_HWI-EAS418:3:7:315:351_length=50#0/1	chr1_14551_14601_+	ENSG00000227232.4
# SRR064286.5709675_HWI-EAS418:3:7:315:351_length=50#0/2	chr1_14681_14731_-	ENSG00000227232.4
# SRR064286.7043546_HWI-EAS418:3:53:1606:39_length=50#0/1	chr1_14912_14962_+	unannotated
# SRR064286.7043546_HWI-EAS418:3:53:1606:39_length=50#0/2 chr1_14928_14978_-	ENSG00000227232.4

# Usage example:
################
# awk -f select_overlappingGene_contiguousMapping.awk <( gzip -dc $outDir/contiguousAlignments_intersected.txt.gz ) 

{	 	
	# Set the mapping block id as the read name
	readId=$4;
	
	# Read mapping coordinates
	chr=$1;
	beg=($2+1);
	end=$3;
	strand=$6;	

	mapCoord=chr"_"beg"_"end"_"strand;
		
	# A) Read do not overlapping any annotated exon
	if ($NF == "0")
	{
	 	gnId="unannotated"; 
	}
		
	# B) Read overlapping an annotated exon
	else
	{
		# Extract gene id
		split($10,feature,":");
		split(feature[2],ids,";");
		gnId=ids[1];
	};

	# Print output:
	row=readId"\t"mapCoord"\t"gnId;
	
	print row;
}
