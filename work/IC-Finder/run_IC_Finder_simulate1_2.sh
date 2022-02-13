#!/bin/bash
datasets=(simulate1 simulate2) 

noiseList=( 0.04 0.08 0.12 0.16 0.20 ) 

for dataset in ${datasets[@]}; do
mkdir -p DOMAINS/${dataset}
for noise in ${noiseList[@]}; do

nice matlab -nodisplay -r  "inputfile='../simulate_data/${dataset}/sim_ob_${noise}.chr5';outputfile='DOMAINS/${dataset}/';filename='${dataset}_${noise}_chr5';",< run_IC_Finder.m

done
done

