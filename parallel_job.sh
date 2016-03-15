#!/bin/bash

cmnd=$1
argfile=$2
progdir=$3
rsltdir=$4

job_system=bsub
que=l

mkdir -p $rsltdir
rm -rf $rsltdir/*

i=0
while read args;do
    jobname="job"$RANDOM
    ./submit_job.sh $job_system $que $i $jobname "${progdir}/${cmnd} $args" 
#    echo $job_system $que $i $jobname "${progdir}/${cmnd} $args" 
    i=`expr $i + 1`
done < $argfile

./monitor

mv par_* ${rsltdir}/.