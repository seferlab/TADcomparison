#!/bin/bash 

#downsample
BIN=50000
echo "Processing size ${BIN}"

ratios=(1 2 5)

mkdir -p DOMAINS/downsample
mkdir -p DATA/downsample_armatus

for ratio in ${ratios[@]}; do
cp ../GM12878_downsample_diff_reso/50k_KR_downsample_ratio_0.0${ratio}.chr6 DATA/downsample_armatus
gzip DATA/downsample_armatus/50k_KR_downsample_ratio_0.0${ratio}.chr6
armatus-2.2/src/armatus -m -r ${BIN} -c 6 -i DATA/downsample_armatus/50k_KR_downsample_ratio_0.0${ratio}.chr6.gz -g 0.5 -o DOMAINS/downsample/50k_KR_downsample_ratio_0.0${ratio}.chr6
done

