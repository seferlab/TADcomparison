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
mkdir -p DATA/${dataset}_armatus

#noise index array
noiseList=(0.04 0.08 0.12 0.16 0.20) 

#start noise index loop
for noise in ${noiseList[@]}; do
cp ../simulate_data/${dataset}/sim_ob_${noise}.chr5 DATA/${dataset}_armatus
mv DATA/${dataset}_armatus/sim_ob_${noise}.chr5 DATA/${dataset}_armatus/${dataset}_${noise}.chr5

gzip DATA/${dataset}_armatus/${dataset}_${noise}.chr5

armatus-2.2/src/armatus -m -r ${BIN} -c 5 -i ./DATA/${dataset}_armatus/${dataset}_${noise}.chr5.gz -g 0.5 -o DOMAINS/${dataset}/${dataset}.${noise}noise

done

done

