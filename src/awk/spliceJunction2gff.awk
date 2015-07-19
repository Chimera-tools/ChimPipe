#!/usr/bin/env awk

# *****************************************************************************
	
#	spliceJunction2gff.awk
	
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

# Takes as input a file with splice junctions produced with "make_spliceJunctions.sh" and produce for each junction two gff features (one for each junction side).

# input
# chr3_182637125_+:chr3_182637126_+ 1 1 1 0 182637098 182637150 GT AG SRR201779.3607803_PATHBIO-SOLEXA2_30TUEAAXX:3:56:1701:1083_length=53#0/1,
# chr4_71629386_+:chr4_71630192_+ 2 2 2 0 71629354 71630235 GT AG SRR201779.54315_PATHBIO-SOLEXA2_30TUEAAXX:3:1:1271:382_length=53#0/1,SRR201779.1335758_PATHBIO-SOLEXA2_30TUEAAXX:3:20:515:1865_length=53#0/1,
# chr20_6855539_+:chr20_24943838_+ 1 1 0 1 6855536 24943886 GC AG SRR201779.5539015_PATHBIO-SOLEXA2_30TUEAAXX:3:86:392:198_length=53#0/2,

# output (GFF)
# chr3	ChimPipe	spliceJunc1	182637098	182637125	.	+	.	JuncId: "chr3_182637125_+:chr3_182637126_+";
# chr3	ChimPipe	spliceJunc2	182637126	182637150	.	+	.	JuncId: "chr3_182637125_+:chr3_182637126_+";
# chr4	ChimPipe	spliceJunc1	71629354	71629386	.	+	.	JuncId: "chr4_71629386_+:chr4_71630192_+";
# chr4	ChimPipe	spliceJunc2	71630192	71630235	.	+	.	JuncId: "chr4_71629386_+:chr4_71630192_+";
# chr20	ChimPipe	spliceJunc1	6855536	6855539	.	+	.	JuncId: "chr20_6855539_+:chr20_24943838_+";
# chr20	ChimPipe	spliceJunc2	24943838	24943886	.	+	.	JuncId: "chr20_6855539_+:chr20_24943838_+";


# Usage example:
################
# awk -f spliceJunction2gff.awk spliceJunc_nbStag_nbtotal_NbUnique_nbMulti_ss1_ss2_sameChrStr_okGxOrder_dist_readIds.txt

{
	split($1,a,":"); 
	split(a[1],b1,"_"); 
	split(a[2],b2,"_"); 
	
	chrA=b1[1]; 
	coordA1=b1[2]; 
	coordA2=$6; 
	strandA=b1[3]; 
	
	chrB=b2[1]; 
	coordB1=b2[2]; 
	coordB2=$7; 
	strandB=b2[3]; 
	
	### Write junction blocks coordinates in such way the begining is lower than the end
	## Block 1
	if (coordA1 < coordA2) 
	{
		begA=coordA1; 
		endA=coordA2;
	}
	else
	{
		begA=coordA2; 
		endA=coordA1;
	} 
	
	## Block 2
	if (coordB1 < coordB2) 
	{
		begB=coordB1; 
		endB=coordB2;
	}
	else
	{
		begB=coordB2; 
		endB=coordB1;
	}
	
	### Make a gff feature for each splice junction staggered block 
	part1=chrA"\tChimPipe\tspliceJunc1\t"begA"\t"endA"\t.\t"strandA"\t.\tJuncId: \""$1"\"\;"; 
	part2=chrB"\tChimPipe\tspliceJunc2\t"begB"\t"endB"\t.\t"strandB"\t.\tJuncId: \""$1"\"\;"; 

	print part1; 
	print part2;
}

