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

# Description
################ 
# Filters chimeric junction candidates to produce a final set of chimeric junctions. For now, the filter cannot be tunned, it filters out chimeric junctions that:
# - Connect genes with high sequence homology (blast alignment length > 30 and %similarity >80)
# - Junction do not supported by at least 10 staggered reads or by at least 1 staggered and 1 pair of reads 

# We will make filter tunable in future releases


# Usage
#######

# awk -f chimjunc_filter.awk chimeric_junctions_candidates.txt


NR==1{
	print $0;
}
(NR>1)&&($0!~/chrM/){

	nbStaggered=$2;
	PElist=$24;
	maxSim=$25;
	maxLgal=$26;
		
	split(PElist,PE,","); 
	maxPE=0; 
	for (nbPE in PE)
	{
		split(PE[nbPE],nb,":");
		if (nb[2]>maxPE)															# Maximum number of encompassing read pairs connecting exons from two different genes
		{
			maxPE=nb[2];
		}
	} 
	
	if (((nbStaggered>=5)&&(maxPE>=0))&&((maxSim<=80)||(maxLgal<=30)||(maxSim==".")))			# First condition
	{
		print $0;																					
	}
	else if (((nbStaggered>=1)&&(maxPE>=1))&&((maxSim<=80)||(maxLgal<=30)||(maxSim==".")))	# Second possible condition
	{
		print $0;
	}
}
