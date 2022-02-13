#!/bin/bash
datasets=(simulate1 simulate2) 

noiseList=( 0.04 0.08 0.12 0.16 0.20 ) 

for dataset in ${datasets[@]}; do
mkdir -p DOMAINS/${dataset}
for noise in ${noiseList[@]}; do
./OnTAD ../simulate_data/${dataset}/sim_ob_${noise}.chr5 -o DOMAINS/${dataset}/${dataset}_${noise}.chr5
done
done


