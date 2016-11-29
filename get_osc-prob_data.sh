#!/bin/bash

### Define get_data function
##### This function runs run.sh to get oscillation probability data for MH and dCP.
get_data () {
    Krho=$1
    MH=$2
    dCP=$3
    if [ $MH -eq 1 ];then
	run_name=NH_CP$dCP
    elif [ $MH -eq -1 ];then
	run_name=IH_CP$dCP
    else
	echo "ERROR: MH value is invalid. 1 or -1 are expected."
    fi
    ./makedir.sh $dir/$run_name 1
    sed -e "s/  MH  1/  MH  $MH/" \
	-e "s/  dCP      0/  dCP      $dCP/" \
	-e "s/   Krho  0d0/   Krho  $Krho/" \
	-e "s/  fdCP          0/  fdCP          $dCP/" params.card_osc-prob-test > params.card
    ./run.sh $run_name 0 1 0 1
    outdir=rslt_$run_name
    mv $outdir/data/prob_* $dir/$run_name
}

### Prepare output directory
#dir=osc-prob_check_vacuum
dir=osc-prob_check_matter
./makedir.sh $dir 1

# take data
Krho=2.9
MH=1
get_data $Krho $MH 0
get_data $Krho $MH 90
get_data $Krho $MH 180
get_data $Krho $MH 270

MH=-1
get_data $Krho $MH 0
get_data $Krho $MH 90
get_data $Krho $MH 180
get_data $Krho $MH 270