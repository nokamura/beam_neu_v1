#!/bin/bash

OAB=2
MH=1
$command ${place} ${beam_mode} ${OAB} ${Y} ${MH} ${dCP} $disp
mv rslt_${place}_${beam_mode}_${OAB}_${Y}_${MH}_${dCP} ${run_dir}
dchi2_20nh=`cat ${run_dir}/rslt_${place}_${beam_mode}_${OAB}_${Y}_${MH}_${dCP}/data/dchi2.dat`
MH=-1
$command ${place} ${beam_mode} ${OAB} ${Y} ${MH} ${dCP} $disp
mv rslt_${place}_${beam_mode}_${OAB}_${Y}_${MH}_${dCP} ${run_dir}
dchi2_20ih=`cat ${run_dir}/rslt_${place}_${beam_mode}_${OAB}_${Y}_${MH}_${dCP}/data/dchi2.dat`
OAB=2.5
MH=1
$command ${place} ${beam_mode} ${OAB} ${Y} ${MH} ${dCP} $disp
mv rslt_${place}_${beam_mode}_${OAB}_${Y}_${MH}_${dCP} ${run_dir}
dchi2_25nh=`cat ${run_dir}/rslt_${place}_${beam_mode}_${OAB}_${Y}_${MH}_${dCP}/data/dchi2.dat`
MH=-1
$command ${place} ${beam_mode} ${OAB} ${Y} ${MH} ${dCP} $disp
mv rslt_${place}_${beam_mode}_${OAB}_${Y}_${MH}_${dCP} ${run_dir}
dchi2_25ih=`cat ${run_dir}/rslt_${place}_${beam_mode}_${OAB}_${Y}_${MH}_${dCP}/data/dchi2.dat`
OAB=3
MH=1
$command ${place} ${beam_mode} ${OAB} ${Y} ${MH} ${dCP} $disp
mv rslt_${place}_${beam_mode}_${OAB}_${Y}_${MH}_${dCP} ${run_dir}
dchi2_30nh=`cat ${run_dir}/rslt_${place}_${beam_mode}_${OAB}_${Y}_${MH}_${dCP}/data/dchi2.dat`
MH=-1
$command ${place} ${beam_mode} ${OAB} ${Y} ${MH} ${dCP} $disp
mv rslt_${place}_${beam_mode}_${OAB}_${Y}_${MH}_${dCP} ${run_dir}
dchi2_30ih=`cat ${run_dir}/rslt_${place}_${beam_mode}_${OAB}_${Y}_${MH}_${dCP}/data/dchi2.dat`

echo $dCP $dchi2_20nh $dchi2_20ih $dchi2_25nh $dchi2_25ih $dchi2_30nh $dchi2_30ih >> dchi2_${place}.tmp