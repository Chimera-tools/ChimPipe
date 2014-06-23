# ~/Awk/junctions_bio_order.awk

BEGIN {while (getline < fileRef >0)
	{
		part1[$1]=$2;
		part2[$1]=$3;
		if (($2=="GT")||($2=="AG"))
		{
			strand1[$1]="+";
		}
		else
		{
			if (($2=="AC")||($2=="CT"))
			{
				strand1[$1]="-";
			}
			else
			{
				strand1[$1]=".";
			}
		}
		if (($3=="GT")||($3=="AG"))
		{
			strand2[$1]="+";
		}
		else
		{
			if (($3=="AC")||($3=="CT"))
			{
				strand2[$1]="-";
			}
			else
			{
				strand2[$1]=".";
			}
		}
	}
}
#   Donor,   Aceptor
#+  GT-AG or GC-AG 
#-  AC-CT or GC-CT
{
	split($1,a,":");
	split(a[1],a1,"_");
	split(a[2],a2,"_");
	if (((part1[$1]=="AG")||(part1[$1]=="CT")||(part2[$1]=="GT")||(part2[$1]=="AC")||(part2[$1]=="GC")) && (stranded==0))
	{
		junc=a2[1]"_"a2[2]"_"strand2[$1]":"a1[1]"_"a1[2]"_"strand1[$1];
		beg=$5
		end=$4
		ss1=part2[$1];
		ss2=part1[$1];
		gnlist1=$10;
		gnlist2=$9;
		gnname1=$12;
		gnname2=$11;
		bt1=$14;
		bt2=$13;
	
	}
	else
	{
		if (stranded==0)
		{
			junc=a1[1]"_"a1[2]"_"strand1[$1]":"a2[1]"_"a2[2]"_"strand2[$1];
		}
		else
		{
			junc=a1[1]"_"a1[2]"_"a1[3]":"a2[1]"_"a2[2]"_"a2[3];
		}
		beg=$4
		end=$5
		ss1=part1[$1];
		ss2=part2[$1];
		gnlist1=$9;
		gnlist2=$10;
		gnname1=$11;
		gnname2=$12;
		bt1=$13;
		bt2=$14;
	}
	print junc, $2, $3, beg, end, $6, $7, $8, ss1, ss2, gnlist1, gnlist2, gnname1, gnname2, bt1, bt2;
}
