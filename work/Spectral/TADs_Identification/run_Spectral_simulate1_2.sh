#!/bin/bash 

checkMakeDirectory(){
        echo -e "checking directory: $1"
        if [ ! -e "$1" ]; then
                echo -e "\tmakedir $1"
                mkdir -p "$1"
        fi
}
checkMakeDirectory ../DOMAINS

datasets=(simulate1 simulate2)
for dataset in ${datasets[@]}; do
mkdir -p ../DOMAINS/${dataset}

#noise index array
noiseList=(0.04 0.08 0.12 0.16 0.20) 

#start noise index loop
for noise in ${noiseList[@]}; do
matlab -nodisplay -r  "file_in='../../simulate_data/${dataset}/sim_ob_${noise}.chr5';file_out='../DOMAINS/${dataset}/${dataset}.${noise}noise'",<Spectral_run2.m
done
done

