#!/bin/bash

### Define compare_data function
compare_data () {
    dir_base=osc-prob_check_base
    dir_check=osc-prob_check_matter
    dir=$1
    detector=$2
    mode=$3

    echo "Location $detector"
    echo "Checking $mode mode"

    file=prob_${mode}_$detector.dat
    file_base=$dir_base/$dir/$file
    file_check=$dir_check/$dir/$file
    
### Extract the 2nd column of $file_check
    cat $file_check | awk '{print $2}' > file_temp
    
### Combine $file_base and the 2nd column of $file_check
    paste -d' ' $file_base file_temp > file_combine
    rm -rf file_temp
    
### check differences between oscillation probabilities in $file_check and $file_base are zero
    flag_diff=0
    nline=1
    while read line; do
	osc1=`echo $line | cut -d ' ' -f 2`
	osc2=`echo $line | cut -d ' ' -f 3`
#    diff=`echo "scale=5; ($osc1 -$osc2)/$osc1" | bc`
#    X=`echo "scale=5; $diff > 0.00001" | bc`
#    if [ $X -eq 1 ];then
	if [ $osc1 != $osc2 ];then
	    echo "ERROR: Line $nline: $line"	
	    flag_diff=1
	    break
	fi
	nline=`expr $nline + 1`
    done < file_combine
    
    if [ $flag_diff -eq 0 ];then
	echo "$mode mode OK"
    fi
    rm -rf file_combine
}

### wrapper of compare_data
compare_dataset () {
    detector=Kr
    dir=$1
    compare_data $dir $detector nm.ne
    compare_data $dir $detector nm.nm
    compare_data $dir $detector nm.nt
    compare_data $dir $detector ne.ne
    compare_data $dir $detector ne.nm
    compare_data $dir $detector ne.nt
    compare_data $dir $detector am.ae
    compare_data $dir $detector am.am
    compare_data $dir $detector am.at
    compare_data $dir $detector ae.ae
    compare_data $dir $detector ae.am
    compare_data $dir $detector ae.at
}

############## Man Program ###################
### prepare oscillation probability data set
./get_osc-prob_data.sh

### check oscillation probability data set
compare_dataset NH_CP0

