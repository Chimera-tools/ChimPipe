# /users/rg/brodriguez/Chimeras_project/Chimeras_detection_pipeline/Chimera_mapping/Workdir/Awk/similarity_filter.awk
# 17/01/2014 


BEGIN{
	while (getline < fileRef >0)
	{
		sim[$1"-"$2]=$3;
		lgal[$1"-"$2]=$4;
	}
}
NR==1{
	print $0, "sim"
}
NR>=2{
	maxLgal=0;
	maxLgalSim=0;
	split($7,gnlist1,",");
	split($8,gnlist2,",");
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
	
#	print "minimun length : ", minLgal
#	print "length alignment: ", maxLgal
	
#	print "maximun similarity : ", simThld
#	print "similarity: ", maxLgalSim
	 
	if (((maxLgal>=minLgal) && (maxLgalSim>=maxSim)))
	{
		print $0, 1
#		print "FILTERED", $0
	}else
	{
		print $0, 0
	}
}



