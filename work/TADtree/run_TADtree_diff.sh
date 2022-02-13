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

ratios=($(seq 1 9))

for ratio in ${ratios[@]}; do
mkdir -p DOMAINS/downsample_${ratio}
nohup python TADtree.py control/GM12878_downsample_${ratio}_chr6.txt >LOG/downsample_${ratio}_chr6.log 2>&1 &
done

#resolutions=(25000 50000 100000)
#for resolution in ${resolutions[@]}; do
#display_reso=`expr $(($resolution/1000))`
#nohup python TADtree.py control/GM12878_${display_reso}k_chr6.txt >LOG/${display_reso}k_chr6.log 2>&1 &
#done
mkdir -p DOMAINS/Gm12878_25k
mkdir -p DOMAINS/Gm12878_50k
nohup python TADtree.py control/GM12878_50k_chr6.txt >LOG/50k_chr6.log 2>&1 &

mkdir -p DOMAINS/Gm12878_100k
nohup python TADtree.py control/GM12878_100k_chr6.txt >LOG/100k_chr6.log 2>&1 &
