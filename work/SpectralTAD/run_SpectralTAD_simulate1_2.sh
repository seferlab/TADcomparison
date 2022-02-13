#!/bin/bash
datasets=(simulate1 simulate2)

noiseList=( 0.04 0.08 0.12 0.16 0.20 )

dataset=simulate1
mkdir -p DOMAINS/${dataset}
for noise in ${noiseList[@]}; do

/home/LiuKun/R-3.6.0/bin/Rscript SpectralTAD_run_non.r -i ../simulate_data/${dataset}/sim_ob_${noise}.chr5 -o DOMAINS/${dataset}/${dataset}_${noise}.chr5 -r 40000 -c chr5
done

dataset=simulate2
mkdir -p DOMAINS/${dataset}
for noise in ${noiseList[@]}; do

/home/LiuKun/R-3.6.0/bin/Rscript SpectralTAD_run.r -i ../simulate_data/${dataset}/sim_ob_${noise}.chr5 -o DOMAINS/${dataset}/${dataset}_${noise}.chr5 -r 40000 -c chr5
done

