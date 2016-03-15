#!/bin/bash

dirCC=$1
dirNC=$2

CC_flag=0
NC_flag=0
if [ $dirCC != "0" ];then
    rm -rf xsecCC
    cp -rf $dirCC xsecCC
    CC_flag=1
else
    echo "skip to set xsecCC"
fi
if [ $dirNC != "0" ];then
    rm -rf xsecNC
    cp -rf $dirNC xsecNC
    NC_flag=1
else
    echo "skip to set xsecNC"
fi

if [ $CC_flag == 1 -o $NC_flag == 1 ];then
    make clean
    make
fi