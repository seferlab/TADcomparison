#!/bin/bash
datasets=(simulate1 simulate2) 

noiseList=( 0.04 0.08 0.12 0.16 0.20 ) 

for dataset in ${datasets[@]}; do
mkdir -p  DOMAINS/${dataset}
for noise in ${noiseList[@]}; do

rows=$(cat ../simulate_data/${dataset}/sim_ob_${noise}.chr5 | wc -l)
let maxn=rows/10
./CHDF ../simulate_data/${dataset}/sim_ob_${noise}.chr5 DOMAINS/${dataset}/${dataset}_${noise}.chr5 ${rows} ${rows} ${maxn}

done
done



