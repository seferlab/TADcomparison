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
checkMakeDirectory control

#downsample
BIN=50000
echo "Processing size ${BIN}"

ratios=(1 2 5)

mkdir -p DOMAINS/downsample
mkdir -p DATA/downsample_MrTADFinder
mkdir -p control/downsample_MrTADFinder

for ratio in ${ratios[@]}; do
python3 convert_triple_MrTADFinder.py ../GM12878_downsample_diff_reso/50k_KR_downsample_ratio_0.0${ratio}.chr6 DATA/downsample_MrTADFinder/50k_KR_downsample_ratio_0.0${ratio}.chr6
nrows=$(cat ../GM12878_downsample_diff_reso/50k_KR_downsample_ratio_0.0${ratio}.chr6 | wc -l)
python3 generate_coor.py control/downsample_MrTADFinder/ratio_0.0${ratio}.file1 control/downsample_MrTADFinder/ratio_0.0${ratio}.file2 ${BIN} ${nrows} chr6
/home/LiuKun/julia/julia-1.4.2/bin/julia run_MrTADFinder.jl DATA/downsample_MrTADFinder/50k_KR_downsample_ratio_0.0${ratio}.chr6 control/downsample_MrTADFinder/ratio_0.0${ratio}.file1 control/downsample_MrTADFinder/ratio_0.0${ratio}.file2 res=1.5 6 DOMAINS/downsample/50k_KR_downsample_ratio_0.0${ratio}.chr6
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
python3 convert_triple_MrTADFinder.py ../GM12878_downsample_diff_reso/50k_KR_downsample_ratio_${ratio}.chr6 DATA/downsample_MrTADFinder/50k_KR_downsample_ratio_${ratio}.chr6
nrows=$(cat ../GM12878_downsample_diff_reso/50k_KR_downsample_ratio_${ratio}.chr6  | wc -l)
python3 generate_coor.py control/downsample_MrTADFinder/ratio_${ratio}.file1 control/downsample_MrTADFinder/ratio_${ratio}.file2 ${BIN} ${nrows} chr6
/home/LiuKun/julia/julia-1.4.2/bin/julia run_MrTADFinder.jl DATA/downsample_MrTADFinder/50k_KR_downsample_ratio_${ratio}.chr6 control/downsample_MrTADFinder/ratio_${ratio}.file1 control/downsample_MrTADFinder/ratio_${ratio}.file2 res=1.5 6 DOMAINS/downsample/50k_KR_downsample_ratio_${ratio}.chr6

done

#different resolutions
mkdir -p DOMAINS/diff_reso
mkdir -p DATA/diff_reso_MrTADFinder
mkdir -p control/diff_reso_MrTADFinder

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
python3 convert_triple_MrTADFinder.py ../GM12878_downsample_diff_reso/${display_reso}k_KR_total.chr6 DATA/diff_reso_MrTADFinder/${display_reso}k_KR_total.chr6
nrows=$(cat ../GM12878_downsample_diff_reso/${display_reso}k_KR_total.chr6 | wc -l)
python3 generate_coor.py control/diff_reso_MrTADFinder/${display_reso}k.file1 control/diff_reso_MrTADFinder/${display_reso}k.file2 ${resolution} ${nrows} chr6
/home/LiuKun/julia/julia-1.4.2/bin/julia run_MrTADFinder.jl DATA/diff_reso_MrTADFinder/${display_reso}k_KR_total.chr6 control/diff_reso_MrTADFinder/${display_reso}k.file1 control/diff_reso_MrTADFinder/${display_reso}k.file2 res=1.5 6 DOMAINS/diff_reso/${display_reso}k_KR_total.chr6
done
