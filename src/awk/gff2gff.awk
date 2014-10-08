#!/usr/bin/env awk

# *****************************************************************************
	
#	gff2gff.awk
	
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
# Needs to be done..

$1!~/#/{
    for (i=1;i<=7;i++)
    {
	printf $i"\t";
    }
    printf $8;
    if(NF>8)
    {
	printf "\t"$9;
	for (i=10;i<=NF;i++)
	{
	    printf " "$i;
	}
    }
    print "";
}

$1~/#/{print}
