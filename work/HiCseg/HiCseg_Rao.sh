#!/bin/bash 

checkMakeDirectory(){
        echo -e "checking directory: $1"
        if [ ! -e "$1" ]; then
                echo -e "\tmakedir $1"
                mkdir -p "$1"
        fi
}
checkMakeDirectory DOMAINS

runHiCseg(){
BIN=$1
echo "Processing size ${BIN}"

list=($(seq $2 $3))
chromList=($(seq 1 22)) 
chromList[${#chromList[*]}]=X

for li in ${list[@]}; do

data=`printf "HIC%03d" $li`

mkdir -p DOMAINS/${data}

for chrom in ${chromList[@]}; do

Rscript hicseg_run_s.r -i ../Rao/${data}/${data}_50k_KR.chr${chrom} -o DOMAINS/${data}/${data}.${chrom}

done
done
}
runHiCseg 50000 1 29
runHiCseg 50000 50 56
runHiCseg 50000 69 74
