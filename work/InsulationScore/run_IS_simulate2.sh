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
#insulation square
square=2000000
#insulation delta span
span=800000

echo "Processing size ${BIN}"

dataset=simulate2

mkdir -p DOMAINS/${dataset}
mkdir -p DATA/${dataset}_insulation

#noise index array
noiseList=(0.04 0.08 0.12 0.16 0.20) 
#start noise index loop
for noise in ${noiseList[@]}; do

python convert_insulation.py -i ../simulate_data/${dataset}/sim_ob_${noise}.chr5 -o DATA/${dataset}_insulation/${dataset}_${noise}.chr5.matrix -c chr5 -b ${BIN}
perl scripts/matrix2insulation.pl -i DATA/${dataset}_insulation/${dataset}_${noise}.chr5.matrix -is ${square} -ids ${span} -im mean -bmoe 3 -nt 0.1 -v	
mv ${dataset}_${noise}.chr5.* ./DOMAINS/${dataset}/


done

