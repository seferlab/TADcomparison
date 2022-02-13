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

#chromosome index array
chromList=($(seq 1 100)) 
chromList[${#chromList[*]}]=X

#start chromosome index loop
for chrom in ${chromList[@]}; do
/home/LiuKun/R-3.6.0/bin/Rscript CHAC_run.r -i ../benchmarks/maptest_${chrom}.txt -o DOMAINS/${dataset}/${dataset}_CHAC.chr${chrom}

done


