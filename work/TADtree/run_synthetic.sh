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

mkdir -p DOMAINS/simulate1
nohup python TADtree.py control/simulate1_control.txt >LOG/simulate1.log 2>&1 &
mkdir -p DOMAINS/simulate2
nohup python TADtree.py control/simulate2_control.txt >LOG/simulate2.log 2>&1 &
mkdir -p DOMAINS/benchmarks
nohup python TADtree.py control/synthetic_data_1_100_control.txt >LOG/benchmarks.log 2>&1 &

