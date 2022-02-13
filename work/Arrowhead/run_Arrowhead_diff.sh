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

chrom=6
ratios=($(seq 1 9))
mkdir -p DOMAINS/downsample

for ratio in ${ratios[@]}; do
java -jar juicer_tools.jar arrowhead /home/chengzj/downsample/preprocess/GM12878/total_merged_nodups_downsample_ratio_${ratio}.hic DOMAINS/downsample/50k_KR_downsample_ratio_${ratio}.chr6 -c chr${chrom} -r ${BIN} -k KR --ignore_sparsity
done

#different resolutions
mkdir -p DOMAINS/diff_reso

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
java -jar juicer_tools.jar arrowhead /home/chengzj/downsample/preprocess/GM12878/total_merged.hic DOMAINS/diff_reso/${display_reso}k_KR_total.chr6 -c chr${chrom} -r ${resolution} -k KR --ignore_sparsity
done
