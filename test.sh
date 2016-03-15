#!/bin/bash

dir=$(cd $(dirname $0);pwd)
rslt=prjct_par
command=job_dchi2.sh
arg[0]="sk 1 3 5 -1 90 0"
arg[1]="ok 1 3 5 -1 90 0"
nargs=1
argfile=args.dat
touch $argfile

i=0
while [ $i -le $nargs ];do
    if [ $i -eq 0 ]; then
	echo ${arg[$i]} > $argfile
    else
	echo ${arg[$i]} >> $argfile
    fi
    i=`expr $i + 1`
done

./parallel_job.sh $command $argfile $dir $rslt 

rm -rf $argfile