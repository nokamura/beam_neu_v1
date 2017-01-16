#!/bin/bash

function check_diff () {
    oab=$1
    file=flux${nu}_$oab.dat
#    diff=`diff $dir_test/$file $dir_ref/flux${nu}_05.dat`
    diff=`diff $dir_test/$file $dir_ref/$file`
    if [ -n "$diff" ];then
	echo $file "Changed"
	ichange=1
    fi
}

function check_diff_nu () {
check_diff 00
check_diff 05
check_diff 06
check_diff 08
check_diff 09
check_diff 10
check_diff 11
check_diff 12
check_diff 13
check_diff 14
check_diff 15
check_diff 20
check_diff 23
check_diff 25
check_diff 30
}

ichange=0

dir_test=rslt_test/data
dir_ref=rslt_base/data

nu=n
check_diff_nu 
nu=a
check_diff_nu 

if [ $ichange -eq 0 ];then
    echo "All test passed"
fi