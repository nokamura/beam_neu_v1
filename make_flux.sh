#!/bin/bash

neu=$1
OAB=$2
name=$3

outfile=$name.inc

sed -e "s/NNEU/$neu/" \
    -e "s/OOAB/$OAB/" ./beam/make_flux.temp > make_flux.f
make flux
./flux

sed -e "s/OOAB/$OAB/" \
    -e "s/oabout/$name/" oabout.inc > beam/$outfile

rm -rf flux make_flux.f oabout.inc