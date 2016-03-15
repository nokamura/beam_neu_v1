#!/bin/bash

dir=$1
#rm -rf xsecCC
rm -rf xsecNC
#rm -rf xsecCC xsecNC
#cp -rf xsecCC_def xsecCC
cp -rf $dir xsecNC
#cp -rf $dir xsecCC
make clean
make