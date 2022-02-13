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

runEAST(){
list=($(seq $1 $2))
BIN=$3
chromList=($(seq 1 22)) 
chromList[${#chromList[*]}]=X

for li in ${list[@]}; do

data=`printf "HIC%03d" $li`

mkdir -p DOMAINS/${data}
mkdir -p DATA/${data}_CaTCH

for chrom in ${chromList[@]}; do
python EAST2.py ../Rao/${data}/${data}_50k_KR.chr${chrom}  DOMAINS/${data}/${data}_EAST.chr${chrom} ${BIN}

done
done
}

runEAST 1 29 50000
runEAST 50 56 50000
runEAST 69 74 50000
