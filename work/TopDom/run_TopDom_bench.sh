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


dataset=bench
mkdir -p DOMAINS/${dataset}
mkdir -p DATA/${dataset}_topdom

#chromosome index array
chromList=($(seq 1 100)) 
chromList[${#chromList[*]}]=X

#start chromosome index loop
for chrom in ${chromList[@]}; do

python convert_N_N+3.py -i ../benchmarks/maptest_${chrom}.txt -o DATA/${dataset}_topdom/${dataset}_topdom.chr${chrom} -c chr${chrom} -b ${BIN}
Rscript topdom_run.r -i DATA/${dataset}_topdom/${dataset}_topdom.chr${chrom} -o DOMAINS/${dataset}/${dataset}_topdom.chr${chrom}


done

