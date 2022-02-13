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
BIN=40000
echo "Processing size ${BIN}"


datasets=(simulate1 simulate2)
for dataset in ${datasets[@]}; do
mkdir -p DOMAINS/${dataset}
mkdir -p DATA/${dataset}_topdom

#noise index array
noiseList=(0.04 0.08 0.12 0.16 0.20) 
#start noise index loop
for noise in ${noiseList[@]}; do

python convert_N_N+3.py -i ../simulate_data/${dataset}/sim_ob_${noise}.chr5 -o DATA/${dataset}_topdom/${dataset}_${noise}.chr5 -c chr5 -b ${BIN}
Rscript topdom_run.r -i DATA/${dataset}_topdom/${dataset}_${noise}.chr5 -o DOMAINS/${dataset}/${dataset}.${noise}noise

done
done

