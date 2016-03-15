#!/bin/bash
    arg[0]="$place $beam_mode 2 $Y 1 0 0"
    arg[1]="$place $beam_mode 2 $Y 1 90 0"
    arg[2]="$place $beam_mode 2 $Y 1 180 0"
    arg[3]="$place $beam_mode 2 $Y 1 270 0"
    arg[4]="$place $beam_mode 2 $Y -1 0 0"
    arg[5]="$place $beam_mode 2 $Y -1 90 0"
    arg[6]="$place $beam_mode 2 $Y -1 180 0"
    arg[7]="$place $beam_mode 2 $Y -1 270 0"
    arg[8]="$place $beam_mode 2.5 $Y 1 0 0"
    arg[9]="$place $beam_mode 2.5 $Y 1 90 0"
    arg[10]="$place $beam_mode 2.5 $Y 1 180 0"
    arg[11]="$place $beam_mode 2.5 $Y 1 270 0"
    arg[12]="$place $beam_mode 2.5 $Y -1 0 0"
    arg[13]="$place $beam_mode 2.5 $Y -1 90 0"
    arg[14]="$place $beam_mode 2.5 $Y -1 180 0"
    arg[15]="$place $beam_mode 2.5 $Y -1 270 0"
    arg[16]="$place $beam_mode 3 $Y 1 0 0"
    arg[17]="$place $beam_mode 3 $Y 1 90 0"
    arg[18]="$place $beam_mode 3 $Y 1 180 0"
    arg[19]="$place $beam_mode 3 $Y 1 270 0"
    arg[20]="$place $beam_mode 3 $Y -1 0 0"
    arg[21]="$place $beam_mode 3 $Y -1 90 0"
    arg[22]="$place $beam_mode 3 $Y -1 180 0"
    arg[23]="$place $beam_mode 3 $Y -1 270 0"
    narg=${#arg[*]}
    i=0
    while [ $i -lt $narg ];do
	echo ${arg[$i]} >> $argfile
	i=`expr $i + 1`
    done