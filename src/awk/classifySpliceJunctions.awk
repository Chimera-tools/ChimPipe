#!/usr/bin/env awk

# *****************************************************************************
	
#	classifySpliceJunctions.awk
	
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

# Takes as input a file containing annotated splice junctions and select normal junctions or chimeric ones (controled with the chimera parameter). A chimeric junction is defined as a junction connecting two different genes

# input
# chr3_42701320_+:chr3_42702977_+	5	6	4	2	42701276	42703014	GT	AG	1	1	1657	100	100	0	0	ENSG00000114853.9	ENSG00000114853.9	ZBTB47	ZBTB47	protein_coding	protein_coding	SRR201779.2086482_PATHBIO-SOLEXA2_30TUEAAXX:3:31:970:1103_length=53#0/1,SRR201779.5570454_PATHBIO-SOLEXA2_30TUEAAXX:3:86:423:861_length=53#0/2,SRR201779.5570454_PATHBIO-SOLEXA2_30TUEAAXX:3:86:423:861_length=53#0/2,SRR201779.5058273_PATHBIO-SOLEXA2_30TUEAAXX:3:79:152:958_length=53#0/2,SRR201779.652134_PATHBIO-SOLEXA2_30TUEAAXX:3:10:1049:1322_length=53#0/2,SRR201779.1924659_PATHBIO-SOLEXA2_30TUEAAXX:3:29:1578:386_length=53#0/2,
# chr1_234492542_+:chr19_1440368_-	4	8	0	8	234492508	1440346	GT	AG	0	na	na	100	100	345	125	ENSG00000235605.1	ENSG00000115268.5	RP5-827C21.1	RPS15	pseudogene	protein_coding	SRR201779.6536907_PATHBIO-SOLEXA2_30TUEAAXX:3:100:478:412_length=53#0/1,SRR201779.5801553_PATHBIO-SOLEXA2_30TUEAAXX:3:89:198:484_length=53#0/1,SRR201779.5773450_PATHBIO-SOLEXA2_30TUEAAXX:3:89:199:486_length=53#0/1,SRR201779.5707857_PATHBIO-SOLEXA2_30TUEAAXX:3:88:814:1476_length=53#0/1,SRR201779.664462_PATHBIO-SOLEXA2_30TUEAAXX:3:10:517:1231_length=53#0/1,SRR201779.2759605_PATHBIO-SOLEXA2_30TUEAAXX:3:42:819:1276_length=53#0/1,SRR201779.3355183_PATHBIO-SOLEXA2_30TUEAAXX:3:52:1539:1691_length=53#0/2,SRR201779.3625018_PATHBIO-SOLEXA2_30TUEAAXX:3:56:1181:711_length=53#0/1,

### output
## normal junctions if "chimera=0"
# chr3_42701320_+:chr3_42702977_+ NORMAL	5	6	4	2	42701276	42703014	GT	AG	1	1	1657	100	100	0	0	ENSG00000114853.9	ENSG00000114853.9	ZBTB47	ZBTB47	protein_coding	protein_coding	SRR201779.2086482_PATHBIO-SOLEXA2_30TUEAAXX:3:31:970:1103_length=53#0/1,SRR201779.5570454_PATHBIO-SOLEXA2_30TUEAAXX:3:86:423:861_length=53#0/2,SRR201779.5570454_PATHBIO-SOLEXA2_30TUEAAXX:3:86:423:861_length=53#0/2,SRR201779.5058273_PATHBIO-SOLEXA2_30TUEAAXX:3:79:152:958_length=53#0/2,SRR201779.652134_PATHBIO-SOLEXA2_30TUEAAXX:3:10:1049:1322_length=53#0/2,SRR201779.1924659_PATHBIO-SOLEXA2_30TUEAAXX:3:29:1578:386_length=53#0/2,

# chr1_234492542_+:chr19_1440368_- CHIMERIC	4	8	0	8	234492508	1440346	GT	AG	0	na	na	100	100	345	125	ENSG00000235605.1	ENSG00000115268.5	RP5-827C21.1	RPS15	pseudogene	protein_coding	SRR201779.6536907_PATHBIO-SOLEXA2_30TUEAAXX:3:100:478:412_length=53#0/1,SRR201779.5801553_PATHBIO-SOLEXA2_30TUEAAXX:3:89:198:484_length=53#0/1,SRR201779.5773450_PATHBIO-SOLEXA2_30TUEAAXX:3:89:199:486_length=53#0/1,SRR201779.5707857_PATHBIO-SOLEXA2_30TUEAAXX:3:88:814:1476_length=53#0/1,SRR201779.664462_PATHBIO-SOLEXA2_30TUEAAXX:3:10:517:1231_length=53#0/1,SRR201779.2759605_PATHBIO-SOLEXA2_30TUEAAXX:3:42:819:1276_length=53#0/1,SRR201779.3355183_PATHBIO-SOLEXA2_30TUEAAXX:3:52:1539:1691_length=53#0/2,SRR201779.3625018_PATHBIO-SOLEXA2_30TUEAAXX:3:56:1181:711_length=53#0/1,


# Usage example:
################
# awk -f classifySpliceJunctions.awk input.txt

BEGIN{
	header="juncCoord\ttype\tnbTotal\tnbStag\tpercStag\tnbMulti\tpercMulti\toverlapA\toverlapB\tdistExonBoundaryA\tdistExonBoundaryB\tdonorSS\tacceptorSS\tbeg\tend\tsameChrStr\tokGxOrder\tdist\tgnIdA\t gnIdB\tgnNameA\tgnNameB\tgnTypeA\tgnTypeB\treadIds";
	print header;
}
{		
	juncCoord=$1;
	nbTotal=$2;
	nbStag=$3;
	percStag=$4;
	nbMulti=$5;
	percMulti=$6;
	overlapA=$7;
	overlapB=$8;
	distExonBoundaryA=$9;
	distExonBoundaryB=$10;
	donorSS=$11;
	acceptorSS=$12;
	beg=$13;
	end=$14;
	sameChrStr=$15;
	okGxOrder=$16;
	dist=$17;
	gnIdA=$18;
	gnIdB=$19;	
	gnNameA=$20;
	gnNameB=$21;
	gnTypeA=$22;
	gnTypeB=$23;
	readIds=$24;
	
	# Set chimeric splice junction as default
	type="CHIMERIC"; 
	
	# A) At least one of junction blocks do not overlap any annotated exon
	if ((gnIdA == "unannotated") || (gnIdB == "unannotated"))
	{
		type="UNANNOTATED"; 
	}
	# B) Both junction blocks overlapping annotated exons
	else
	{
		# If gene name not available use gene id
		if ((gnNameA == "na") || (gnNameB == "na"))
		{
			nbA=split(gnIdA, gnListA, ","); 
			nbB=split(gnIdB,gnListB, ",");
		}
		# If gene name available use gene name (better, since for gencode there are several different gene ids for some gene within the same annotation)
		else 
 		{
			nbA=split(gnNameA, gnListA, ","); 
			nbB=split(gnNameB,gnListB, ",");
		}
	
		# Do all the possible gene pairs comparisons between the 2 list of genes overlapping junction donor and acceptor respectivelly. 
		# If the donor and acceptor overlap the same gene, set type to normal junction. Otherwise, leave it as chimeric junction (default)
		for (i=1;i<=nbA;i++)
		{
			for (j=1;j<=nbB;j++)
			{
				if (gnListA[i]==gnListB[j])
				{
					type="NORMAL"; 
				}
			}
		}
	}
		
	## Report splice junctions + type + associated info:
	row=juncCoord"\t"type"\t"nbTotal"\t"nbStag"\t"percStag"\t"nbMulti"\t"percMulti"\t"overlapA"\t"overlapB"\t"distExonBoundaryA"\t"distExonBoundaryB"\t"donorSS"\t"acceptorSS"\t"beg"\t"end"\t"sameChrStr"\t"okGxOrder"\t"dist"\t"gnIdA"\t" gnIdB"\t"gnNameA"\t"gnNameB"\t"gnTypeA"\t"gnTypeB"\t"readIds;
	print row;  
}


