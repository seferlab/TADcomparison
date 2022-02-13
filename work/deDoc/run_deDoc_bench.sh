#!/bin/bash

chromList=($(seq 1 100)) 

dataset=bench
mkdir -p DOMAINS/${dataset}
for chrom in ${chromList[@]}; do
python3 convert_triple.py ../benchmarks/maptest_${chrom}.txt DOMAINS/${dataset}/${dataset}_${chrom}.txt
java -jar deDoc.jar  DOMAINS/${dataset}/${dataset}_${chrom}.txt
done

