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
/home/LiuKun/R-3.6.2/bin/Rscript TADBD_run.r -i ../benchmarks/maptest_${chrom}.txt -o DOMAINS/${dataset}/${dataset}_TADBD.chr${chrom}

done


