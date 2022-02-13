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

ratios=($(seq 1 9))

mkdir -p DOMAINS/downsample

for ratio in ${ratios[@]}; do
/home/LiuKun/R-3.6.0/bin/Rscript CHAC_run.r -i ../GM12878_downsample_diff_reso/50k_KR_downsample_ratio_${ratio}.chr6 -o DOMAINS/downsample/50k_KR_downsample_ratio_${ratio}.chr6
done

#different resolutions
mkdir -p DOMAINS/diff_reso

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
/home/LiuKun/R-3.6.0/bin/Rscript CHAC_run.r -i ../GM12878_downsample_diff_reso/${display_reso}k_KR_total.chr6 -o DOMAINS/diff_reso/${display_reso}k_KR_total.chr6
done
