#!/usr/bin/env awk

# Description
###############
# Takes as input CRAC output file for chimeric junction (--chimera option has to be enabled) and produces chimeric junctions in ChimPipe's format plus their support (nb. reads spanning the junction point)

## Input: CRAC output file for chimeric junctions
#read_id tag_chimera loc_end_first_exon_on_genome loc_start_second_exon_on_genome pos_junction_on_read single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc
# 29456 duplicate chr4|1,76807216 chr2|-1,191402746 pos_junction=23 chr4|1,76807192 pos_location=3 GACAGGTTAGTTTTACCCTACTGATGATGTGTTGTTGCCATGGTAATCCT 20318,20391,21933,21663,21662,25106,24953,24585,

## Output: Chimeric junctions in ChimPipe's format + number of reads supporting each junction
# chimJunc nbSupporingReads
# chr4_76807216_+:chr2_191402746_-

# Usage:
########
# awk -f make_chimJunjc_crac_output.awk crac_output_chimeras


BEGIN{
	print "# chimJunc", "nbSupporingReads";
}
! /^#/{
	### 5 prime junction boundary
	junc5prime=$3;
	gsub(/\|/, ",", junc5prime);
	split(junc5prime, a,",");
	chr5prime=a[1];
	coord5prime=a[3];
	
	# Convert CRAC strand information (1=="+"; -1=="-") into +,-
	strand5prime=(a[2]=="1" ? "+" : "-");
	
	### 3 prime junction boundary
	junc3prime=$4;
	gsub(/\|/, ",", junc3prime);
	split(junc3prime, b,",");
	chr3prime=b[1];
	coord3prime=b[3];
	
	# Convert CRAC strand information (1=="+"; -1=="-") into +,-
	strand3prime=(b[2]=="1" ? "+" : "-");
	
	### Make chimJunc
	chimJunc=chr5prime"_"coord5prime"_"strand5prime":"chr3prime"_"coord3prime"_"strand3prime

	supportSpanningReads[chimJunc]=supportSpanningReads[chimJunc]+1	
}

END{
	for (junc in supportSpanningReads)
	{
		nbSpanningReads=supportSpanningReads[junc];
		print junc, nbSpanningReads;
	}
}
