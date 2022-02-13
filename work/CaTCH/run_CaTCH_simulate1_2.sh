#!/bin/bash 

checkMakeDirectory(){
        echo -e "checking directory: $1"
        if [ ! -e "$1" ]; then
                echo -e "\tmakedir $1"
                mkdir -p "$1"
        fi
}
checkMakeDirectory DOMAINS
checkMakeDirectory DATA

datasets=(simulate1 simulate2)
for dataset in ${datasets[@]}; do
mkdir -p DOMAINS/${dataset}
mkdir -p DATA/${dataset}_CaTCH

#noise index array
noiseList=(0.04 0.08 0.12 0.16 0.20) 

#start noise index loop
for noise in ${noiseList[@]}; do
python3 convert_triple_CaTCH.py ../simulate_data/${dataset}/sim_ob_${noise}.chr5 DATA/${dataset}_CaTCH/sim_ob_${noise}.chr5 5
gzip DATA/${dataset}_CaTCH/sim_ob_${noise}.chr5

Rscript synthetic_CaTCH_run.r -i DATA/${dataset}_CaTCH/sim_ob_${noise}.chr5.gz -o DOMAINS/${dataset}/${dataset}.${noise}noise

done

done

