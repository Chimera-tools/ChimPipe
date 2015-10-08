#!/usr/bin/env awk

# *****************************************************************************
	
#	cutgff.awk
	
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

# Cuts a gff file to the toth field (outpuf: from field no 1 to no to
# which is specified as an argument)

# awk -v to=10 -f ~/Awk/cutgff.awk in.gff > out.gff


{
    s="";
    if(to<=9)
    {
	for(i=1; i<=to-1; i++)
	{
	    s=(s)($i)("\t");
	}  
	s=(s)($i);
	print s;
    }
    else
    {
	for(i=1; i<=8; i++)
	{
	    s=(s)($i)("\t");
	}  
	for(i=9; i<=to-1; i++) 
	{
	    s=(s)($i)" ";
	}	
	s=(s)($i);
	print s;
    }
}
