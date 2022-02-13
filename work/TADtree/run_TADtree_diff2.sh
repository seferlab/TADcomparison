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

ratios=(1 2 5)

for ratio in ${ratios[@]}; do
mkdir -p DOMAINS/downsample_0.0${ratio}
nohup python TADtree.py control/GM12878_downsample_0.0${ratio}_chr6.txt >LOG/downsample_0.0${ratio}_chr6.log 2>&1 &
done

