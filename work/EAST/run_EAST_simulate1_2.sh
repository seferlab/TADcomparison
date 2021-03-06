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

#noise index array
noiseList=(0.04 0.08 0.12 0.16 0.20) 

#start noise index loop
for noise in ${noiseList[@]}; do
python3 EAST2.py ../simulate_data/${dataset}/sim_ob_${noise}.chr5  DOMAINS/${dataset}/${dataset}.${noise}noise  40000

done

done

