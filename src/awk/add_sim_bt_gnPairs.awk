#!/usr/bin/env awk

# *****************************************************************************
	
#	add_sim_bt_gnPairs.awk
	
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
# Need to be done..

BEGIN{
	while (getline < fileRef >0)
	{
		sim[$1"-"$2]=$3;
		lgal[$1"-"$2]=$4;
	}
}

# Print header
NR==1{
	print $0, "maxSim", "maxLgal";
}

NR>1{
	maxLgal=0;
	maxLgalSim=0;
	split($13,gnlist1,",");
	split($14,gnlist2,",");
	for (gn1 in gnlist1)
	{
		for (gn2 in gnlist2)
		{
			if((gnlist1[gn1]!="")||(gnlist2[gn2]!=""))
			{
				gp=((gnlist1[gn1]<=gnlist2[gn2]) ? (gnlist1[gn1]"-"gnlist2[gn2]) : (gnlist2[gn2]"-"gnlist1[gn1]));
				if (lgal[gp]>maxLgal)
				{
					maxLgal=lgal[gp];
					maxLgalSim=sim[gp];
				}
			}
		}
	}
	if((maxLgalSim!=0)||(maxLgal!=0))	 
	{
		print $0, maxLgalSim, maxLgal; 
	}
	else
	{
		print $0, ".", ".";
	} 
}



