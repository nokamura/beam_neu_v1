#!/bin/bash

function check_diff () {
    nu1=$1
    nu2=$2
    file=prob_$nu1.${nu2}_${detector}.dat
#    diff=`diff $dir_test/$file $dir_ref/flux${nu}_05.dat`
    diff=`diff $dir_test/$file $dir_ref/$file`
    if [ -n "$diff" ];then
	echo $file "Changed"
	echo $diff
	ichange=1
    fi
}

function check_diff_nu () {
nu_type=$1
nu1=$2

nu=$nu_type$nu1
check_diff $nu ${nu_type}e
check_diff $nu ${nu_type}m
check_diff $nu ${nu_type}t
}

#########################################################################
####### MAIN PROGRAM
#########################################################################
ichange=0

dir_test=rslt_test/data
dir_ref=rslt_base/data

detector=SK
check_diff_nu n e
check_diff_nu n m
check_diff_nu a e
check_diff_nu a m

if [ $ichange -eq 0 ];then
    echo "All test passed"
fi