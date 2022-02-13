#!/bin/bash 

checkMakeDirectory(){
        echo -e "checking directory: $1"
        if [ ! -e "$1" ]; then
                echo -e "\tmakedir $1"
                mkdir -p "$1"
        fi
}
checkMakeDirectory ../DOMAINS

dataset=bench
mkdir -p ../DOMAINS/${dataset}

chromList=($(seq 1 100)) 

#start chromosome index loop
for chrom in ${chromList[@]}; do
matlab -nodisplay -r  "file_in='../../benchmarks/maptest_${chrom}.txt';file_out='../DOMAINS/${dataset}/${dataset}_Spectral.chr${chrom}'",<Spectral_run.m
done


