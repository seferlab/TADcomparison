#!/bin/bash
datasets=(simulate1 simulate2 simulate3 simulate4) 

noiseList=( 0.04 0.08 0.12 0.16 0.20 ) 

for dataset in ${datasets[@]}; do
mkdir -p ../DOMAINS/${dataset}
for noise in ${noiseList[@]}; do

nice matlab -nodisplay -r  "filepath='../../simulate_data/${dataset}/';name='sim_ob_${noise}.chr5';Res=40000;chromo='${noise}';outputfolder_name='../DOMAINS/${dataset}';Max_TADsize='600000';",< ClusterTAD_main.m

done
done

