#!/bin/bash

dataset=bench

chromList=($(seq 1 100)) 
#chromList[${#chromList[*]}]=X

mkdir -p DOMAINS/${dataset}
for chrom in ${chromList[@]}; do

rows=$(cat ../benchmarks/maptest_${chrom}.txt | wc -l)
let maxn=rows/10
./CHDF ../benchmarks/maptest_${chrom}.txt DOMAINS/${dataset}/${dataset}.${chrom} ${rows} ${rows} ${maxn}

done


