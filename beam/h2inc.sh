#!/bin/bash

file=$1
name=${file%.h}

while read line;do
    echo C$line >> temp1
done < $file

i=1
while read line;do
    head=${line:0:2}
    if [ $head = "C{" ];then
	echo $i$line >> temp2
	i=`expr $i + 1`
    else
	echo $line >> temp2
    fi
done < temp1

sed -e "s%\(.*\)C{\(.*\)},%      data (${name}(\1,oi),oi=1,5) / \2 /%" \
    -e "s%\(.*\)C{\(.*\)}%      data (${name}(\1,oi),oi=1,5) / \2 /%" \
    -e "/Cd/d" \
    -e "/C}/d" temp2 > temp3

sed -e "/-\*/a GG1" temp3 > temp4
sed -e "s/GG1/      real*8 ${name}(300,5)/" temp4 > ${name}.inc

rm -rf temp*