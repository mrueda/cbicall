#!/bin/bash
set -eu 

dir=/media/mrueda/2TBS/CNAG/Project_CBI_Call
cbicall=/media/mrueda/2TBS/CNAG/Project_CBI_Call/cbicall/bin/cbicall
ncpu=4

for dirname in MA99999_exome
do
 cd $dir/$dirname
 echo $dirname
 for sample in MA*{P,M,F}*ex
 do
  echo "...$sample"
  cd $sample
  cat<<EOF>$sample.smk_wes_single.yaml
mode: single
pipeline: wes
workflow_engine: snakemake
sample: $dir/$dirname/$sample
EOF
time $cbicall -t $ncpu -p $sample.smk_wes_single.yaml > $sample.smk_wes_single.log 2>&1
  cd ..
 done
cd ..
done
