#!/bin/bash

dataset=bench

chromList=($(seq 1 100)) 
#chromList[${#chromList[*]}]=X

mkdir -p DOMAINS/${dataset}
for chrom in ${chromList[@]}; do

nice matlab -nodisplay -r  "inputfile='../benchmarks/maptest_${chrom}.txt';outputfile='DOMAINS/${dataset}/';filename='${dataset}_${chrom}';",< run_IC_Finder.m

done




