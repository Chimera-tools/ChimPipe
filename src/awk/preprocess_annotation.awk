#!/usr/bin/env awk

# *****************************************************************************
	
#	preprocess_annotation.awk
	
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
# Takes as input a raw gencode gene annotation file in gff format and report a gff with a non-reduntant list of exons (one per row) with basic information associated (genomic coordinates, gene id, gene name and gene type)  

### input
# chr1	HAVANA	exon	11869	12227	.	+	.	gene_id "ENSG00000223972.4"; transcript_id "ENST00000456328.2"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_status "KNOWN"; transcript_name "DDX11L1-002"; exon_number 1; exon_id "ENSE00002234944.1"; level 2; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
# chr1	HAVANA	exon	12613	12721	.	+	.	gene_id "ENSG00000223972.4"; transcript_id "ENST00000456328.2"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_status "KNOWN"; transcript_name "DDX11L1-002"; exon_number 2; exon_id "ENSE00003582793.1"; level 2; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";

### output
# chr1	HAVANA	exon	11869	12227	.	+	.	gene_id "ENSG00000223972.4"; gene_name "DDX11L1";  gene_type "pseudogene"; 
# chr1	ENSEMBL	exon	12613	12721	.	+	.	gene_id "ENSG00000223972.4"; gene_name "DDX11L1"; gene_type "pseudogene"; 


# Usage example:
################
# awk  -f preprocess_annotation.awk gencode19_annot.gff

$3=="exon"{
	exon=$1"_"$4"_"$5"_"$7;  
	source[exon]=$2; 
	score[exon]=$6; 
	frame[exon]=$8; 
	addGene="1";
	
	for(i=9; i<=(NF); i++)
	{	
		if ($i=="gene_id")
		{
			split($(i+1), a,"\"" );
			
			# A) No gene id associated to the exon -> Asign gene id to the exon
			if (gene_id[exon]=="")
			{
				gene_id[exon]=a[2];
			}
			# B) Already at least one gene id associated to the exon
			else
			{
				# Check if the gene id is already in the list of the genes associated to the exon
				split(gene_id[exon], geneIdList, ",");
				for (id in geneIdList)
				{
					# Gene id already in the list -> do not add it
					if (geneIdList[id] == a[2])
					{
						addGene="0";
					}
				}
				
				# Gene not in the list -> add it
				if (addGene=="1")
				{
					gene_id[exon]=gene_id[exon]","a[2];
				}
			}
			break
		}
	}
	
	for(i=9; i<=(NF); i++)
	{		
		## A) Add gene type if gene not in the list of genes associated to the exon
		if (($i=="gene_type")&&(addGene=="1"))
		{
			split($(i+1), a,"\"" );
			
			if (gene_type[exon]=="")
			{
				gene_type[exon]=a[2];
			}
			else
			{
				gene_type[exon]=gene_type[exon]","a[2];
			}
		}
		## Add gene name if gene not in the list of genes associated to the exon
		else if (($i=="gene_name")&&(addGene=="1"))
		{			
			split($(i+1), a,"\"" );
			
			if (gene_name[exon]=="")
			{
				gene_name[exon]=a[2];
			}
			else
			{
				gene_name[exon]=gene_name[exon]","a[2];
			}
		}
	}
}

END{
	for (exon in source)
	{
		split(exon,pos,"_"); 
		chr=pos[1]; 
		beg=pos[2]; 
		end=pos[3]; 
		strand=pos[4]; 
		
		# A) If gene_name information not provided use "na"
		if (gene_name[exon]=="")
		{
			gene_name[exon]="na"; 
		}
		
		# B) If gene_type information not provided use "na"
		if (gene_type[exon]=="")
		{
			gene_type[exon]="na";
		}
		
		# Print exon information in gtf format as output
		printf chr"\t"source[exon]"\texon\t"beg"\t"end"\t"score[exon]"\t"strand"\t"frame[exon]"\tgene_id ""\""gene_id[exon]"\"; gene_name ""\""gene_name[exon]"\";"" gene_type ""\""gene_type[exon]"\";\n"; 
	}
}


