#!/bin/bash

dataset=bench

chromList=($(seq 1 100)) 

mkdir ../DOMAINS/${dataset}
for chrom in ${chromList[@]}; do
nice matlab -nodisplay -r  "filepath='../../benchmarks/';name='maptest_${chrom}.txt';Res=40000;chromo='bench${chrom}';outputfolder_name='../DOMAINS/${dataset}';Max_TADsize='600000';",< ClusterTAD_main.m
done

