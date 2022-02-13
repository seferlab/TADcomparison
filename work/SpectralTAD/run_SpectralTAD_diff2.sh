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

for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
/home/LiuKun/R-3.6.0/bin/Rscript SpectralTAD_run.r -i ../GM12878_downsample_diff_reso/50k_KR_downsample_ratio_${ratio}.chr6 -o DOMAINS/downsample/50k_KR_downsample_ratio_${ratio}.chr6 -r ${BIN} -c chr6
done

