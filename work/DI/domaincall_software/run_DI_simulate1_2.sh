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

BIN=40000
echo "Processing size ${BIN}"

datasets=(simulate1 simulate2)
for dataset in ${datasets[@]}; do

mkdir -p DATA/${dataset}_DI
mkdir -p HMM/${dataset}

#chromosome index array
noiseList=(0.04 0.08 0.12 0.16 0.20) 

#start chromosome index loop
for noise in ${noiseList[@]} 
do
fai_path=./special/${dataset}_${noise}.chr5

#convert raw observed symmetric matrix to N*(N+3) matrix
python ../convert_N_N+3.py -i ../../simulate_data/${dataset}/sim_ob_${noise}.chr5 -o DATA/${dataset}_DI/${dataset}_DI.chr5.${noise}noise -c chr5 -b ${BIN}

## Compute Directionality Index for every chromosome
echo "Computing Directionality Index chr${chrom}"
## one .matrix file for each chromosome. contains the binned genome 
perl ${perl_path}/DI_from_matrix.pl DATA/${dataset}_DI/${dataset}_DI.chr5.${noise}noise ${BIN} 2000000 ${fai_path}> DI/${dataset}_chr5.${noise}noise.DI

## Run HMM 
echo "Running HMM for noise${noise}"
file1="DI/${dataset}_chr5.${noise}noise.DI"
file2="DI/${dataset}_chr5.${noise}noise.out"
nice matlab -nodisplay -r  "inputfile='${file1}',outputfile='${file2}'",< HMM_calls_rep.m &> LOG/hmm_${dataset}.${noise}noise.log

## postprocessing:: Clean up and generate 7 col
echo "Post processing files"
perl ${perl_path}/file_ends_cleaner.pl DI/${dataset}_chr5.${noise}noise.out DI/${dataset}_chr5.${noise}noise.DI| perl ${perl_path}/converter_7col.pl > HMM/hmm_${dataset}.${noise}noise.hmm


echo "noise${noise}"
echo "Post processing files: Generating Domains"
perl ${perl_path}/hmm_probablity_correcter.pl HMM/hmm_${dataset}.${noise}noise.hmm 2 0.99 ${BIN} | perl ${perl_path}/hmm-state_caller.pl ${fai_path} chr5 | perl ${perl_path}/hmm-state_domains.pl > DOMAINS/${dataset}.${noise}noise.domain


done


done

