#!/bin/bash 

#downsample
BIN=50000
echo "Processing size ${BIN}"

ratios=(1 2 5)

mkdir -p DOMAINS/downsample

for ratio in ${ratios[@]}; do
/home/LiuKun/R-3.6.0/bin/Rscript CHAC_run.r -i ../GM12878_downsample_diff_reso/50k_KR_downsample_ratio_0.0${ratio}.chr6 -o DOMAINS/downsample/50k_KR_downsample_ratio_0.0${ratio}.chr6
done

