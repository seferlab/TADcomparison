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
max_size=750000

ratios=(1 2 5)

for ratio in ${ratios[@]}; do
mkdir -p ../DOMAINS/downsample_${ratio}
matlab -nodisplay -r  "filepath='../../GM12878_downsample_diff_reso/';name='50k_KR_downsample_ratio_0.0${ratio}.chr6';Res=${BIN};chromo='chr6';outputfolder_name='../DOMAINS/downsample_0.0${ratio}';Max_TADsize=${max_size};",< ClusterTAD_main.m
done

