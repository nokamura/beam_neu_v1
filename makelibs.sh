#!/bin/bash

rm -rf lib
mkdir lib

cd minuit_f
make clean
make
cd ..

cd mylib
make clean
make
cd ..

echo `pwd` > maindir.txt

./set_xsec.sh xsecNC_def