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

#downsample
BIN=50000
echo "Processing size ${BIN}"

ratios=(1 2 5)
mkdir -p DOMAINS/downsample
mkdir -p DATA/downsample

for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
python3 convert_N_N+3.py -i ../GM12878_downsample_diff_reso/50k_KR_downsample_ratio_${ratio}.chr6 -o DATA/downsample/downsample_ratio_${ratio}.chr6 -c chr6 -b ${BIN}
Rscript topdom_run.r -i DATA/downsample/downsample_ratio_${ratio}.chr6 -o DOMAINS/downsample/50k_KR_downsample_ratio_${ratio}.chr6
done
