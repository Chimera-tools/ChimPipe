#!/bin/bash
file=$1

if [[ $file ]]; then
    command="cat"
    if [[ $file =~ .*\.gz ]];then
        command="zcat"
    fi
    command="$command $file | "
fi

command="${command}awk 'BEGIN{for(i=1;i<=256;i++){ord[sprintf(\"%c\",i)]=i;}}NR%4==0{split(\$0,a,\"\");for(i=1;i<=length(a);i++){if(ord[a[i]]<59){print \"Offset 33\";exit 0}if(ord[a[i]]>74){print \"Offset 64\";exit 0}}}'"

eval $command
