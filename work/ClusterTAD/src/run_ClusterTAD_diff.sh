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

ratios=($(seq 1 9))

for ratio in ${ratios[@]}; do
mkdir -p ../DOMAINS/downsample_${ratio}
matlab -nodisplay -r  "filepath='../../GM12878_downsample_diff_reso/';name='50k_KR_downsample_ratio_${ratio}.chr6';Res=${BIN};chromo='chr6';outputfolder_name='../DOMAINS/downsample_${ratio}';Max_TADsize=${max_size};",< ClusterTAD_main.m
done

#different resolutions
mkdir -p ../DOMAINS/diff_reso

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
max_size=`expr $(($resolution*15))`
matlab -nodisplay -r  "filepath='../../GM12878_downsample_diff_reso/';name='${display_reso}k_KR_total.chr6';Res=${resolution};chromo='chr6';outputfolder_name='../DOMAINS/diff_reso';Max_TADsize=${max_size};",< ClusterTAD_main.m
done
