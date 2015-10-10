#!/usr/bin/env awk

# *****************************************************************************
	
#	bamHeader2genomeFile.awk
	
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
##############

# Produce a genome file from a BAM file header

## Usage:
# samtools view -H file.bam | awk -f bamHeader2genomeFile.awk  

## Input: bam header
# @HD	VN:1.3	SO:coordinate
# @SQ	SN:chr1	LN:197195432
# @SQ	SN:chr10	LN:129993255
# @SQ	SN:chr11	LN:121843856
# @SQ	SN:chr12	LN:121257530
# @SQ	SN:chr13	LN:120284312
# @SQ	SN:chr14	LN:125194864
# @SQ	SN:chr15	LN:103494974
# @SQ	SN:chr16	LN:98319150
# @SQ	SN:chr17	LN:95272651
# @SQ	SN:chr18	LN:90772031
# @SQ	SN:chr19	LN:61342430
# @SQ	SN:chr2	LN:181748087
# @SQ	SN:chr3	LN:159599783
# @SQ	SN:chr4	LN:155630120
# @SQ	SN:chr5	LN:152537259
# @SQ	SN:chr6	LN:149517037
# @SQ	SN:chr7	LN:152524553
# @SQ	SN:chr8	LN:131738871
# @SQ	SN:chr9	LN:124076172
# @SQ	SN:chrM	LN:16299
# @SQ	SN:chrX	LN:166650296
# @SQ	SN:chrY	LN:15902555
# @RG	ID:0	PG:GEM	PL:ILLUMINA	SM:0
# @PG	ID:GEM	PN:gem-2-sam	VN:1.847

## Output: genome file
# chr1	197195432
# chr10	129993255
# chr11	121843856
# chr12	121257530
# chr13	120284312
# chr14	125194864
# chr15	103494974
# chr16	98319150
# chr17	95272651
# chr18	90772031
# chr19	61342430
# chr2	181748087
# chr3	159599783
# chr4	155630120
# chr5	152537259
# chr6	149517037
# chr7	152524553
# chr8	131738871
# chr9	124076172
# chrM	16299
# chrX	166650296
# chrY	15902555

# Select @SQ lines
$1 == "@SQ"{
	split($2, seqNameList, ":");
	split($3, seqLengthList, ":");
	
	seqName = seqNameList[2];
	seqLength = seqLengthList[2];
	
	
	row=seqName"\t"seqLength;
	print row;
}	

