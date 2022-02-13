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
mkdir -p DATA/${dataset}

#noise index array
noiseList=(0.04 0.08 0.12 0.16 0.20) 

#start noise index loop
for noise in ${noiseList[@]}; do
python convert_homer.py -i ../simulate_data/${dataset}/sim_ob_${noise}.chr5 -o DATA/${dataset}/sim_ob_${noise}.chr5 -c chr5 -b 40000
/home/chengzj/mypython/bin/hicConvertFormat -m DATA/${dataset}/sim_ob_${noise}.chr5 --inputFormat homer --outputFormat h5 -o DATA/${dataset}/sim_ob_${noise}.chr5.h5
/home/chengzj/mypython/bin/hicFindTADs  -m DATA/${dataset}/sim_ob_${noise}.chr5.h5 --outPrefix DOMAINS/${dataset}/${dataset}.${noise}noise --numberOfProcessors 16 --correctForMultipleTesting None --minDepth 120000 --maxDepth 240000 --step 40000
done
done

