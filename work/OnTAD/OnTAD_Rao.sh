#!/bin/bash 

checkMakeDirectory(){
        echo -e "checking directory: $1"
        if [ ! -e "$1" ]; then
                echo -e "\tmakedir $1"
                mkdir -p "$1"
        fi
}
checkMakeDirectory DOMAINS

runOnTAD(){
BIN=$1
echo "Processing size ${BIN}"

list=($(seq $2 $3))
chromList=($(seq 1 22)) 
chromList[${#chromList[*]}]=X

for li in ${list[@]}; do

data=`printf "HIC%03d" $li`

mkdir -p DOMAINS/${data}

for chrom in ${chromList[@]}; do

./OnTAD ../Rao/${data}/${data}_50k_KR.chr${chrom} -o DOMAINS/${data}/${data}.chr${chrom}

done
done
}
runOnTAD 50000 1 29
runOnTAD 50000 50 56
runOnTAD 50000 69 74
