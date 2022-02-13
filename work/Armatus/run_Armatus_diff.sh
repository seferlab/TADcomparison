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
mkdir -p DATA/downsample_armatus

for ratio in ${ratios[@]}; do
cp ../GM12878_downsample_diff_reso/50k_KR_downsample_ratio_${ratio}.chr6 DATA/downsample_armatus
gzip DATA/downsample_armatus/50k_KR_downsample_ratio_${ratio}.chr6
armatus-2.2/src/armatus -m -r ${BIN} -c 6 -i DATA/downsample_armatus/50k_KR_downsample_ratio_${ratio}.chr6.gz -g 0.5 -o DOMAINS/downsample/50k_KR_downsample_ratio_${ratio}.chr6
done

cp ../GM12878_downsample_diff_reso/50k_KR_total.chr6 DATA/downsample_armatus
gzip DATA/downsample_armatus/50k_KR_total.chr6
armatus-2.2/src/armatus -m -r ${BIN} -c 6 -i DATA/downsample_armatus/50k_KR_total.chr6.gz -g 0.5 -o DOMAINS/downsample/50k_KR_total.chr6

#different resolutions
mkdir -p DOMAINS/diff_reso
mkdir -p DATA/diff_reso_armatus

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
cp ../GM12878_downsample_diff_reso/${display_reso}k_KR_total.chr6 DATA/diff_reso_armatus
gzip DATA/diff_reso_armatus/${display_reso}k_KR_total.chr6
armatus-2.2/src/armatus -m -r ${resolution} -c 6 -i DATA/diff_reso_armatus/${display_reso}k_KR_total.chr6.gz -g 0.5 -o DOMAINS/diff_reso/${display_reso}k_KR_total.chr6
done
