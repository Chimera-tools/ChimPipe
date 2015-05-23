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

# Takes as input a file 

# input

# output

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
		if ($i=="gene_name")
		{
			split($(i+1), a,"\"" );
			
			if (gene_name[exon]=="")
			{
				gene_name[exon]=a[2];
			}
			else
			{
				split(gene_name[exon], geneNameList, ",");
				for (name in geneNameList)
				{
					if (geneNameList[name] == a[2])
					{
						addGene="0";
					}
				}
				
				if (addGene=="1")
				{
					gene_name[exon]=gene_name[exon]","a[2];
				}
			}
			break
		}
	}
	
	for(i=9; i<=(NF); i++)
	{		
		##
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
		else if (($i=="gene_id")&&(addGene=="1"))
		{			
			split($(i+1), a,"\"" );
			
			if (gene_id[exon]=="")
			{
				gene_id[exon]=a[2];
			}
			else
			{
				gene_id[exon]=gene_id[exon]","a[2];
			}
		}
		else if (($i=="transcript_id")&&(addGene=="1"))
		{			
			split($(i+1), a,"\"" );
			
			if (transcript_id[exon]=="")
			{
				transcript_id[exon]=a[2];
			}
			else
			{
				transcript_id[exon]=transcript_id[exon]","a[2];
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
		
		if (gene_id[exon]=="")
		{
			gene_id[exon]="NA"; 
		}
		
		if (gene_type[exon]=="")
		{
			gene_type[exon]="NA";
		}
			
		printf chr"\t"source[exon]"\texon\t"beg"\t"end"\t"score[exon]"\t"strand"\t"frame[exon]"\tgene_name ""\""gene_name[exon]"\";"" gene_id ""\""gene_id[exon]"\";"" gene_type ""\""gene_type[exon]"\";"" transcript_id ""\""transcript_id[exon]"\";\n"; 
	}
}


