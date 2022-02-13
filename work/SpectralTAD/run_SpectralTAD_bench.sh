#!/bin/bash

chromList=($(seq 1 100)) 

dataset=bench
mkdir -p DOMAINS/${dataset}
for chrom in ${chromList[@]}; do
/home/LiuKun/R-3.6.0/bin/Rscript SpectralTAD_run.r -i ../benchmarks/maptest_${chrom}.txt -o DOMAINS/${dataset}/${dataset}_${chrom}.txt -r 40000 -c bench${chrom}
done

