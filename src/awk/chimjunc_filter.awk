#!/usr/bin/env awk

# *****************************************************************************
	
#	chimjunc_filter.awk
	
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

# Usage
#######

# awk -v filterConf="string with filter configuration" chimeric_junctions.txt

# Example: 

# awk -v filterConf="10,2,90,25"


## Input 
########

# 1) raw chimeric junctions file with header
############################################
# Example:
# juncId nbstag nbtotal maxbeg maxEnd samechr samestr dist ss1 ss2 gnlist1 gnlist2 gnname1 gnname2 bt1 bt2 PEsupport maxSim maxLgal
# chr11_60781486_+:chr11_120100337_+ 1 1 60781475 120100375 1 1 59318851 GT AG CD6, OAF, CD6, OAF, . . . . .

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


NR==1{print}
(NR>=2)&&($0!~/chrM/){
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
