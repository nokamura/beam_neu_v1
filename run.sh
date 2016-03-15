#!/bin/bash
if [[ "$1" == "-h" ]]; then
    echo ""
    echo "Usage: run.sh [run_name] [run mode] [display] [mail] [imake]"
    echo ""
    exit
fi
selfdir=$(cd $(dirname $0);pwd)

run=test
run_mode=0
display=1
mail=0
imake=1
if [[ $1 = "h" ]]; then
    echo "input run name"
    read run
    echo "input run mode: 0:Serial 1:Parallel"
    read run_mode
    echo "display"
    read display  
    echo "input whether inform when the job finishes or not 0:no 1:inform"
    read mail
    echo "do make before running? 0:no 1:yes"
    read imake
else
    if [ $# -ge 1 ];then
	run=$1
    fi
    if [ $# -ge 2 ];then
	run_mode=$2
    fi
    if [ $# -ge 3 ];then
	display=$3
    fi
    if [ $# -ge 4 ];then
	mail=$4
    fi
    if [ $# -ge 5 ];then
	imake=$5
    fi
fi    

run_dir=rslt_${run}
outfile=${run_dir}/results.log
defout1=output1
results=results_main
defout2=output2
touch $defout1
touch $defout2
disp=0
if [ $display -eq 1 ]; then
    disp=1
fi
cluster=kekcc

mkdir -p plots
mkdir -p data
rm -rf plots/*
rm -rf data/* 

start_time=`date '+%s'`
date=`date '+%Y/%m/%d'`
ttime=`date '+%T'`
echo ${date} ${ttime} >> $defout1  
echo "" >> $defout1  

if [ -e ${run_dir} ]; then
    rm -rf ${run_dir}/*
    echo "run directory exists." >> $defout1  
    echo " " >> $defout1  
else
    mkdir -p ${run_dir}
fi

### start program ###
if [ $imake -eq 1 ];then
    make
fi
if [ ${run_mode} -eq 0 ]; then  # serial
    ${selfdir}/main < params.card > $results
elif [ ${run_mode} -eq 1 ]; then  # parallel
    echo "Sorry, the parallel mode is not implemented yet." 
    echo "Exiting program..." 
    exit
fi    

mv *.dat data/.
cp -rf data ${run_dir}

### end program ###
    
end_time=`date '+%s'`
SS=`expr ${end_time} - ${start_time}` 
HH=`expr ${SS} / 3600` 
SS=`expr ${SS} % 3600` 
MM=`expr ${SS} / 60` 
SS=`expr ${SS} % 60` 
elapsed_time="${HH}:${MM}:${SS}" 
echo "total time = " $elapsed_time >> ${defout2}
echo "" >> $defout2
echo `date '+%T'` >> $defout2
echo " " >> ${defout2}

if [ $mail -eq 1 ]; then
    if [ $cluster == "kekcc" ]; then
	bsub -q e -J $run -u takaesu@post.kek.jp nulljob.sh >/dev/null 2>&1
    else
	echo "Notification mail cannot be send from this cluster system. Exting..." >> ${defout2}
	echo " " >> ${defout2}
	exit
    fi
fi

cat $defout1 > $outfile 
cat $results >> $outfile
cat $defout2 >> $outfile
if [ $disp -eq 1 ];then
    cat $outfile
    echo " "
fi
cp -rf params.card ${run_dir}/.

rm -rf $defout1 $results $defout2