#!/usr/bin/env awk

# Description
###############
# Takes as input PRADA candidate file for chimeric junction and produces chimeric junctions in ChimPipe's format ((donchr_doncoord_donstrand:accchr_acccoord_accstrand))

# Usage:
########
# awk -f make_chimJunc_prada.awk prada_output_chimeras

### Input: PRADA output file for chimeric junctions (fusim.fus.candidates.txt file)
# Gene_A        Gene_B  A_chr   B_chr   A_strand        B_strand        Discordant_n    JSR_n   perfectJSR_n    Junc_n  Position_Consist        Junction
# AMHR2   TMTC2   12      12      1       1       14      24      24      1       PARTIALLY       AMHR2:12:53818254_TMTC2:12:83358803,24
# GLCCI1  WEE2    7       7       1       1       8       15      15      1       PARTIALLY       GLCCI1:7:8099878_WEE2:7:141420735,15
# TP53INP2        CPNE9   20      3       1       1       5       14      14      1       YES     TP53INP2:20:33297328_CPNE9:3:9754225,14
# 7 (11 fields)
# 198 (12 fields)  

### output: Chimeric junctions in ChimPipe's format
# junc_id
# chr12_53818254_+:chr12_83358803_+
# chr7_8099878_+:chr7_141420735_+


NR==1{
    print "junc_id";
}

((NR>=2)&&(NF==12)){
    split($12,a,"_"); 
    split(a[1],a1,":"); 
    split(a[2],a2,":"); 
    split(a2[3],b,","); 
    print "chr"a1[2]"_"a1[3]"_"($5>0 ? "+" :"-")":chr"a2[2]"_"b[1]"_"($6>0 ? "+" :"-");
}
