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

dataset=bench
mkdir -p DOMAINS/${dataset}

chromList=($(seq 1 100)) 

#start chromosome index loop
for chrom in ${chromList[@]}; do
python3 EAST2.py  ../benchmarks/maptest_${chrom}.txt  DOMAINS/${dataset}/${dataset}_EAST.chr${chrom} 40000

done


