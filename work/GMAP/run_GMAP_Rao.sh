#!/bin/bash 

checkMakeDirectory(){
        echo -e "checking directory: $1"
        if [ ! -e "$1" ]; then
                echo -e "\tmakedir $1"
                mkdir -p "$1"
        fi
}
checkMakeDirectory DOMAINS

runGMAP(){
list=($(seq $1 $2))
BIN=$3
chromList=($(seq 1 22)) 
chromList[${#chromList[*]}]=X

for li in ${list[@]}; do

data=`printf "HIC%03d" $li`

mkdir -p DOMAINS/${data}

for chrom in ${chromList[@]}; do
Rscript GMAP_run.r -i ../Rao/${data}/${data}_50k_KR.chr${chrom} -o DOMAINS/${data}/${data}_GMAP.chr${chrom} -r ${BIN}
done
done
}

runGMAP 1 29 50000
runGMAP 50 56 50000
runGMAP 69 74 50000
