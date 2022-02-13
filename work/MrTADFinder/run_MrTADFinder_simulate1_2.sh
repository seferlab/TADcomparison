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
checkMakeDirectory control

BIN=40000
echo "Processing size ${BIN}"

datasets=(simulate1 simulate2)
for dataset in ${datasets[@]}; do
mkdir -p DOMAINS/${dataset}
mkdir -p DATA/${dataset}_MrTADFinder
mkdir -p control/${dataset}_MrTADFinder

#noise index array
noiseList=(0.04 0.08 0.12 0.16 0.20) 

#start noise index loop
for noise in ${noiseList[@]}; do
python3 convert_triple_MrTADFinder.py ../simulate_data/${dataset}/sim_ob_${noise}.chr5 DATA/${dataset}_MrTADFinder/sim_ob_${noise}.chr5
nrows=$(cat ../simulate_data/${dataset}/sim_ob_${noise}.chr5 | wc -l)
python3 generate_coor.py control/${dataset}_MrTADFinder/${noise}.file1 control/${dataset}_MrTADFinder/${noise}.file2 ${BIN} ${nrows} chr5
/home/LiuKun/julia/julia-1.4.2/bin/julia run_MrTADFinder.jl DATA/${dataset}_MrTADFinder/sim_ob_${noise}.chr5 control/${dataset}_MrTADFinder/${noise}.file1 control/${dataset}_MrTADFinder/${noise}.file2 res=3.0 5 DOMAINS/${dataset}/${dataset}.${noise}noise
done

done

