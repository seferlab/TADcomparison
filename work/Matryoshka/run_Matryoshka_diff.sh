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
mkdir -p DATA/downsample_Matryoshka

for ratio in ${ratios[@]}; do
python3 convert_N_N+3.py -i  ../GM12878_downsample_diff_reso/50k_KR_downsample_ratio_0.0${ratio}.chr6 -o DATA/downsample_Matryoshka/50k_KR_downsample_ratio_0.0${ratio}.chr6 -c chr6 -b 50000
gzip DATA/downsample_Matryoshka/50k_KR_downsample_ratio_0.0${ratio}.chr6
/home/LiuKun/tmp/matryoshka/matryoshka-master/build/src/matryoshka -r ${BIN} -c 6 -i DATA/downsample_Matryoshka/50k_KR_downsample_ratio_0.0${ratio}.chr6.gz -g 0.5 -o DOMAINS/downsample/50k_KR_downsample_ratio_0.0${ratio}.chr6
done

ratios=($(seq 1 9))

for ratio in ${ratios[@]}; do
python3 convert_N_N+3.py -i ../GM12878_downsample_diff_reso/50k_KR_downsample_ratio_${ratio}.chr6 -o DATA/downsample_Matryoshka/50k_KR_downsample_ratio_${ratio}.chr6 -c chr6 -b 50000
gzip DATA/downsample_Matryoshka/50k_KR_downsample_ratio_${ratio}.chr6
/home/LiuKun/tmp/matryoshka/matryoshka-master/build/src/matryoshka -r ${BIN} -c 6 -i DATA/downsample_Matryoshka/50k_KR_downsample_ratio_${ratio}.chr6.gz -g 0.5 -o DOMAINS/downsample/50k_KR_downsample_ratio_${ratio}.chr6
done

#different resolutions
mkdir -p DOMAINS/diff_reso
mkdir -p DATA/diff_reso_Matryoshka

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
python3 convert_N_N+3.py -i ../GM12878_downsample_diff_reso/${display_reso}k_KR_total.chr6 -o DATA/diff_reso_Matryoshka/${display_reso}k_KR_total.chr6 -c chr6 -b ${resolution}
gzip DATA/diff_reso_Matryoshka/${display_reso}k_KR_total.chr6
/home/LiuKun/tmp/matryoshka/matryoshka-master/build/src/matryoshka -r ${resolution} -c 6 -i DATA/diff_reso_Matryoshka/${display_reso}k_KR_total.chr6.gz -g 0.5 -o DOMAINS/diff_reso/${display_reso}k_KR_total.chr6
done
