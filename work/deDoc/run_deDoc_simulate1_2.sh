#!/bin/bash
datasets=(simulate1 simulate2) 

noiseList=( 0.04 0.08 0.12 0.16 0.20) 

for dataset in ${datasets[@]}; do
mkdir -p DOMAINS/${dataset}
for noise in ${noiseList[@]}; do
python3 convert_triple.py ../simulate_data/${dataset}/sim_ob_${noise}.chr5 DOMAINS/${dataset}/${dataset}_${noise}.chr5
java -jar deDoc.jar  DOMAINS/${dataset}/${dataset}_${noise}.chr5
done
done

