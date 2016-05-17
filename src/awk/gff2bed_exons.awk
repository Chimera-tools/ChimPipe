#!/usr/bin/env awk

# *****************************************************************************
	
#	gff2bed_exons.awk
	
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

# Takes as input a version 2 gff file containing annotated exons and report all the information in bed format. Keep gene id, name and biotype information as bed feature name. 

### input:
# chr1	HAVANA	exon	11869	12227	.	+	.	gene_id "ENSG00000223972.4"; gene_name "DDX11L1"; gene_type "pseudogene"; 
# chr1	ENSEMBL	exon	12613	12721	.	+	.	gene_id "ENSG00000223972.4"; gene_name "DDX11L1"; gene_type "pseudogene"; 

### outut:
# chr1	11868	12227	ENSG00000223972.4:exon	.	+
# chr1	12612	12721	ENSG00000223972.4:exon	.	+

# Usage example:
################
# awk  -f gff2bed_exons.awk annotatedExons.[gff2|gtf] 

{
	chr=$1;
	beg=($4-1);
	end=$5;
	feature=$3;
	strand=$7;

	## Parse attribute fields to extract gene_id, gene_name and gene_type information
	for(i=9; i<=(NF); i++)
	{	
		split($(i+1), a,"\"" );

		if ($i=="gene_id")
		{
			gnId=a[2];
		}
		else if ($i=="gene_name")
                {
                        gnName=a[2];
                }
		else if	($i=="gene_type")
                {
                        gnType=a[2];
                }
	}
	
	## Make bed entry name:
	bedEntryName=feature":"gnId";"gnName";"gnType;

	# Make and print output row:
	row=chr"\t"beg"\t"end"\t"bedEntryName"\t.\t"strand;
		
	print row
}


