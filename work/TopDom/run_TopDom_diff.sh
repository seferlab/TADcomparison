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
mkdir -p DATA/downsample

for ratio in ${ratios[@]}; do

python3 convert_N_N+3.py -i ../GM12878_downsample_diff_reso/50k_KR_downsample_ratio_${ratio}.chr6 -o DATA/downsample/downsample_ratio_${ratio}.chr6 -c chr6 -b ${BIN}
Rscript topdom_run.r -i DATA/downsample/downsample_ratio_${ratio}.chr6 -o DOMAINS/downsample/50k_KR_downsample_ratio_${ratio}.chr6
done

#different resolutions
mkdir -p DOMAINS/diff_reso
mkdir -p DATA/diff_reso

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`

python3 convert_N_N+3.py -i ../GM12878_downsample_diff_reso/${display_reso}k_KR_total.chr6 -o DATA/diff_reso/${display_reso}k_KR_total.chr6 -c chr6 -b ${resolution}
Rscript topdom_run.r -i DATA/diff_reso/${display_reso}k_KR_total.chr6 -o DOMAINS/diff_reso/${display_reso}k_KR_total.chr6
done
