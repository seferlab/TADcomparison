#!/bin/bash

chromList=($(seq 1 100)) 

dataset=bench
mkdir -p DOMAINS/${dataset}
for chrom in ${chromList[@]}; do
./OnTAD ../benchmarks/maptest_${chrom}.txt -o DOMAINS/${dataset}/${dataset}_${chrom}.txt
done

