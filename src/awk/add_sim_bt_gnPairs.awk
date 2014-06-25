# /users/rg/brodriguez/Chimeras_project/Chimeras_detection_pipeline/Chimera_mapping/Workdir/Awk/similarity_filter.awk
# 17/01/2014 

BEGIN{
	while (getline < fileRef >0)
	{
		sim[$1"-"$2]=$3;
		lgal[$1"-"$2]=$4;
	}
}
{
	maxLgal=0;
	maxLgalSim=0;
	split($11,gnlist1,",");
	split($12,gnlist2,",");
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



