#!/bin/bash 

checkMakeDirectory(){
        echo -e "checking directory: $1"
        if [ ! -e "$1" ]; then
                echo -e "\tmakedir $1"
                mkdir -p "$1"
        fi
}
checkMakeDirectory DOMAINS
BIN=40000
echo "Processing size ${BIN}"


dataset=bench
mkdir -p DOMAINS/${dataset}

#chromosome index array
chromList=($(seq 1 100)) 

#start chromosome index loop
for chrom in ${chromList[@]}; do
Rscript hicseg_run_s.r -i ../benchmarks/maptest_${chrom}.txt -o DOMAINS/${dataset}/${dataset}_hicseg.chr${chrom}

done

