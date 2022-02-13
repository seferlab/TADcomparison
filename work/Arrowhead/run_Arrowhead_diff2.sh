#!/bin/bash 

#downsample
BIN=50000
echo "Processing size ${BIN}"

ratios=(1 2 5)

mkdir -p DOMAINS/downsample

for ratio in ${ratios[@]}; do
java -jar /root/liukun/downsample/juicer_tools.jar arrowhead /root/liukun/downsample/preprocess/GM12878/total_merged_nodups_downsample_ratio_0.0${ratio}.hic DOMAINS/downsample/50k_KR_downsample_ratio_0.0${ratio}.chr6 -c chr6 -r ${BIN} -k KR --ignore_sparsity
done

