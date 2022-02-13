#!/bin/bash 

checkMakeDirectory(){
        echo -e "checking directory: $1"
        if [ ! -e "$1" ]; then
                echo -e "\tmakedir $1"
                mkdir -p "$1"
        fi
}
checkMakeDirectory DOMAINS
#downsample
ratios=(1 2 5)

mkdir -p DOMAINS/downsample
mkdir -p DATA/downsample

for ratio in ${ratios[@]}; do
python convert_homer.py -i ../GM12878_downsample_diff_reso/50k_KR_downsample_ratio_0.0${ratio}.chr6 -o DATA/downsample/50k_KR_downsample_ratio_0.0${ratio}.chr6 -c chr6 -b 50000
/home/chengzj/mypython/bin/hicConvertFormat -m DATA/downsample/50k_KR_downsample_ratio_0.0${ratio}.chr6 --inputFormat homer --outputFormat h5 -o DATA/downsample/50k_KR_downsample_ratio_0.0${ratio}.chr6.h5
/home/chengzj/mypython/bin/hicFindTADs  -m DATA/downsample/50k_KR_downsample_ratio_0.0${ratio}.chr6.h5 --outPrefix DOMAINS/downsample/50k_KR_downsample_ratio_0.0${ratio}.chr6 --numberOfProcessors 16 --correctForMultipleTesting None --minDepth 150000 --maxDepth 300000 --step 50000
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
python convert_homer.py -i ../GM12878_downsample_diff_reso/50k_KR_downsample_ratio_${ratio}.chr6 -o DATA/downsample/50k_KR_downsample_ratio_${ratio}.chr6 -c chr6 -b 50000
/home/chengzj/mypython/bin/hicConvertFormat -m DATA/downsample/50k_KR_downsample_ratio_${ratio}.chr6 --inputFormat homer --outputFormat h5 -o DATA/downsample/50k_KR_downsample_ratio_${ratio}.chr6.h5
/home/chengzj/mypython/bin/hicFindTADs  -m DATA/downsample/50k_KR_downsample_ratio_${ratio}.chr6.h5 --outPrefix DOMAINS/downsample/50k_KR_downsample_ratio_${ratio}.chr6 --numberOfProcessors 16 --correctForMultipleTesting None --minDepth 150000 --maxDepth 300000 --step 50000
done

#different resolutions
mkdir -p DOMAINS/diff_reso
mkdir -p DATA/diff_reso

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
para1=`expr $(($resolution*3))`
para2=`expr $(($resolution*6))`
python convert_homer.py -i ../GM12878_downsample_diff_reso/${display_reso}k_KR_total.chr6 -o DATA/diff_reso/${display_reso}k_KR_total.chr6 -c chr6 -b ${resolution}
/home/chengzj/mypython/bin/hicConvertFormat -m DATA/diff_reso/${display_reso}k_KR_total.chr6 --inputFormat homer --outputFormat h5 -o DATA/diff_reso/${display_reso}k_KR_total.chr6.h5
/home/chengzj/mypython/bin/hicFindTADs  -m DATA/diff_reso/${display_reso}k_KR_total.chr6.h5 --outPrefix DOMAINS/diff_reso/${display_reso}k_KR_total.chr6 --numberOfProcessors 16 --correctForMultipleTesting None --minDepth ${para1} --maxDepth ${para2} --step ${resolution}
done
