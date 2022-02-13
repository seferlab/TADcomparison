#!/bin/bash 

checkMakeDirectory(){
        echo -e "checking directory: $1"
        if [ ! -e "$1" ]; then
                echo -e "\tmakedir $1"
                mkdir -p "$1"
        fi
}
checkMakeDirectory DOMAINS

dataset=bench
mkdir -p DOMAINS/${dataset}

chromList=($(seq 1 100)) 

#start chromosome index loop
for chrom in ${chromList[@]}; do
Rscript GMAP_run.r -i ../benchmarks/maptest_${chrom}.txt -o DOMAINS/${dataset}/${dataset}_GMAP.chr${chrom} -r 40000

done


