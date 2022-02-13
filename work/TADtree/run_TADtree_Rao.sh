#!/bin/bash

checkMakeDirectory(){
        echo -e "checking directory: $1"
        if [ ! -e "$1" ]; then
                echo -e "\tmakedir $1"
                mkdir -p "$1"
        fi
}
checkMakeDirectory DOMAINS
checkMakeDirectory LOG

runTADtree(){
list=($(seq $1 $2))

for li in ${list[@]}; do

data=`printf "HIC%03d" $li`

mkdir -p DOMAINS/${data}
nohup /home/chengzj/miniconda2/bin/python TADtree.py control/${data}.txt >LOG/${data}.log 2>&1 &
done
}
runTADtree 1 29
runTADtree 50 56
runTADtree 69 74
