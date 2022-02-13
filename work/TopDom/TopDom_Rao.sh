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

runTopDom(){
BIN=$1
echo "Processing size ${BIN}"

list=($(seq $2 $3))
chromList=($(seq 1 22)) 
chromList[${#chromList[*]}]=X

for li in ${list[@]}; do

data=`printf "HIC%03d" $li`

mkdir -p DOMAINS/${data}
mkdir -p DATA/${data}_TopDom

for chrom in ${chromList[@]}; do

python3 convert_N_N+3.py -i ../Rao/${data}/${data}_50k_KR.chr${chrom} -o DATA/${data}_TopDom/${data}_TopDom.chr${chrom} -c chr${chrom} -b ${BIN}
Rscript topdom_run.r -i DATA/${data}_TopDom/${data}_TopDom.chr${chrom} -o DOMAINS/${data}/${data}.chr${chrom}

done
done
}
runTopDom 50000 1 29
runTopDom 50000 50 56
runTopDom 50000 69 74
