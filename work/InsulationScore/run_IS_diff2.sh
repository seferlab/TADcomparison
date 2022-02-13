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
square=2500000
span=1000000
echo "Processing size ${BIN}"

ratios=(1 2 5)
mkdir -p DOMAINS/downsample
mkdir -p DATA/downsample

for ratio in ${ratios[@]}; do
python convert_insulation.py -i ../GM12878_downsample_diff_reso/50k_KR_downsample_ratio_0.0${ratio}.chr6 -o DATA/downsample/downsample_ratio_0.0${ratio}_chr6.matrix -c chr6 -b ${BIN}
perl scripts/matrix2insulation.pl -i DATA/downsample/downsample_ratio_0.0${ratio}_chr6.matrix -is ${square} -ids ${span} -im mean -bmoe 3 -nt 0.1	
mv downsample_ratio_0.0${ratio}_chr6.* DOMAINS/downsample
done

