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
ratios=(1 2 5)

mkdir -p DOMAINS/downsample
mkdir -p DATA/downsample_CaTCH

for ratio in ${ratios[@]}; do
python3 convert_triple_CaTCH.py ../GM12878_downsample_diff_reso/50k_KR_downsample_ratio_0.0${ratio}.chr6 DATA/downsample_CaTCH/50k_KR_downsample_ratio_0.0${ratio}.chr6 6
gzip  DATA/downsample_CaTCH/50k_KR_downsample_ratio_0.0${ratio}.chr6
Rscript CaTCH_run.r -i DATA/downsample_CaTCH/50k_KR_downsample_ratio_0.0${ratio}.chr6.gz -o DOMAINS/downsample/50k_KR_downsample_ratio_0.0${ratio}.chr6
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
python3 convert_triple_CaTCH.py ../GM12878_downsample_diff_reso/50k_KR_downsample_ratio_${ratio}.chr6 DATA/downsample_CaTCH/50k_KR_downsample_ratio_${ratio}.chr6 6
gzip DATA/downsample_CaTCH/50k_KR_downsample_ratio_${ratio}.chr6
Rscript CaTCH_run.r -i DATA/downsample_CaTCH/50k_KR_downsample_ratio_${ratio}.chr6.gz -o DOMAINS/downsample/50k_KR_downsample_ratio_${ratio}.chr6
done

#different resolutions
mkdir -p DOMAINS/diff_reso
mkdir -p DATA/diff_reso_CaTCH

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
python3 convert_triple_CaTCH.py ../GM12878_downsample_diff_reso/${display_reso}k_KR_total.chr6 DATA/diff_reso_CaTCH/${display_reso}k_KR_total.chr6 6
gzip DATA/diff_reso_CaTCH/${display_reso}k_KR_total.chr6
Rscript CaTCH_run.r -i DATA/diff_reso_CaTCH/${display_reso}k_KR_total.chr6.gz -o DOMAINS/diff_reso/${display_reso}k_KR_total.chr6
done
