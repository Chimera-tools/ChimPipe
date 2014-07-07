
# Usage
#######

# awk -v filterConf="string with filter configuration" chimeric_junctions.txt

# Example: 

# awk -v filterConf="10,2,90,25"


## Input 
########

# 1) raw chimeric junctions file
################################
# Example:
# chr10_64973501_.:chr20_62152688_- 1 1 64973482 62152656 0 NA NA NA CT JMJD1C, PPDPF, JMJD1C, PPDPF, . . . . .
# chr20_4162515_+:chr6_32083374_. 1 1 4162478 32083360 0 NA NA GT NA SMOX, ATF6B, SMOX, ATF6B, . . . 100.00 17
# chr16_67144142_+:chr16_67700168_. 1 2 67144118 67700141 0 NA NA GT NA C16orf70, ENKD1, C16orf70, ENKD1, . . 1-1:2, . .
# chr16_67144142_+:chr16_67700168_. 1 2 67144118 67700141 0 NA NA GT NA C16orf70, PNLC1, GATC2, EALC, . . 1-2:24, 90.00 25 
 
# 2) string with the configuration for the filter
#################################################
#Example: 1,2,75,50;

# 4 fields where:
# 	1st: minimum number of staggered reads spanning the chimeric junction.
#	2nd: maximum number of paired-end reads encompassing the chimeric junction.
#   3rd: maximum similarity between the connected genes.
# 	4rd: minimum length of the high similar region between the connected genes. 	

# All these conditions have to be fulfilled for a junction to pass the filter. 

# Setting two different conditions:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# It is also possible to set two different conditions and only one of them has to be fulfilled to pass the filtering:

#  Example: 10,0,80,30;1,1,80,30;

# This means that the junction has to be supported by at least 10 staggered reads or by 1 staggered and pair of reads 
# encompassing the junction point. 


NR>1{
	n=split(filterConf,a,";"); 
	split(a[1],conf1,","); 
	split($17,PE,","); 
	maxPE=0; 
	for (nbPE in PE)
	{
		split(PE[nbPE],nb,":");
		if (nb[2]>maxPE)															# Maximum number of encompassing read pairs connecting exons from two different genes
		{
			maxPE=nb[2];
		}
	} 
	minStag1=conf1[1]
	minPE1=conf1[2]; 
	maxSim1=conf1[3];
	minLgAl1=conf1[4];  
	if (($2>=minStag1)&&(maxPE>=minPE1)&&(($18<=maxSim1)||($19<=minLgAl1)))			# First condition
	{
		filtered="0";																# If first condition is fulfilled don't filter out the junction						
	}
	else
	{
		filtered="1";
		if (n=="3")																	# If there are two conditions... 10,0,80,30;1,1,80,30;
		{ 
			split(a[2],conf2,","); 
			minStag2=conf2[1]
			minPE2=conf2[2]; 
			maxSim2=conf2[3]; 
			minLgAl2=conf2[4];
			if (($2>=minStag2)&&(maxPE>=minPE2)&&(($18<=maxSim2)||($19<=minLgAl2))) 	# Second condition
			{
				filtered="0";														# If second condition is fulfilled don't filter out the junction
			}		
		} 	
	}
	if (filtered=="0")
	{
		print $0;
	}
}
