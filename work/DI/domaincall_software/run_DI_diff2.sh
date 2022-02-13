#!/bin/bash

checkMakeDirectory(){
        echo -e "checking directory: $1"
        if [ ! -e "$1" ]; then
                echo -e "\tmakedir $1"
                mkdir -p "$1"
        fi
}
# create subdirectories
checkMakeDirectory DI
checkMakeDirectory HMM
checkMakeDirectory DOMAINS
checkMakeDirectory LOG
checkMakeDirectory DATA

perl_path=./perl_scripts
fai_path=../hg19.fai

BIN=50000
echo "Processing size ${BIN}"

mkdir -p DATA/downsample_DI

ratios=(1 2 5)

for ratio in ${ratios[@]}; do
#convert raw observed symmetric matrix to N*(N+3) matrix
python ../convert_N_N+3.py -i ../../GM12878_downsample_diff_reso/50k_KR_downsample_ratio_0.0${ratio}.chr6 -o DATA/downsample_DI/downsample_0.0${ratio}.chr6 -c chr6 -b ${BIN}

## Compute Directionality Index for every chromosome
## one .matrix file for each chromosome. contains the binned genome 
perl ${perl_path}/DI_from_matrix.pl DATA/downsample_DI/downsample_0.0${ratio}.chr6 ${BIN} 2500000 ${fai_path}> DI/downsample_0.0${ratio}.chr6.DI

## Run HMM 
echo "Running HMM for downsample_0.0${ratio}"
file1="DI/downsample_0.0${ratio}.chr6.DI"
file2="DI/downsample_0.0${ratio}.chr6.out"
nice matlab -nodisplay -r  "inputfile='${file1}',outputfile='${file2}'",< HMM_calls.m &> LOG/hmm_downsample_0.0${ratio}.log

## postprocessing:: Clean up and generate 7 col
echo "Post processing files"
perl ${perl_path}/file_ends_cleaner.pl DI/downsample_0.0${ratio}.chr6.out DI/downsample_0.0${ratio}.chr6.DI| perl ${perl_path}/converter_7col.pl > HMM/hmm_downsample_0.0${ratio}.hmm


echo "downsample_0.0${ratio}"
echo "Post processing files: Generating Domains"
perl ${perl_path}/hmm_probablity_correcter.pl HMM/hmm_downsample_0.0${ratio}.hmm 2 0.99 ${BIN} | perl ${perl_path}/hmm-state_caller.pl ${fai_path} chr6 | perl ${perl_path}/hmm-state_domains.pl > DOMAINS/downsample_0.0${ratio}.chr6.domain

done

