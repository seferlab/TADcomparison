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
fai_path=../bench.fai
BIN=10000
echo "Processing size ${BIN}"



dataset=bench

mkdir -p DATA/${dataset}_DI
mkdir -p HMM/${dataset}

chromList=($(seq 1 100)) 

for chrom in ${chromList[@]} ; do
#convert raw observed symmetric matrix to N*(N+3) matrix
python ../convert_N_N+3.py -i ../../benchmarks/maptest_${chrom}.txt -o DATA/${dataset}_DI/${dataset}_DI.chr${chrom} -c chr${chrom} -b ${BIN}

## Compute Directionality Index for every chromosome
echo "Computing Directionality Index chr${chrom}"
## one .matrix file for each chromosome. contains the binned genome 
perl ${perl_path}/DI_from_matrix.pl DATA/${dataset}_DI/${dataset}_DI.chr${chrom} ${BIN} 100000 ${fai_path}> DI/${dataset}_chr${chrom}.DI

done


## concatenate all the files together

echo "Concatenating all files together. "
cat /dev/null > DI/${dataset}.DI
for chrom in ${chromList[@]}; do
cat DI/${dataset}_chr${chrom}.DI >> DI/${dataset}.DI
done

## Run HMM 
echo "Running HMM for chr${chrom}"
file1="DI/${dataset}.DI"
file2="DI/${dataset}.out"
nice matlab -nodisplay -r  "inputfile='${file1}',outputfile='${file2}'",< HMM_calls.m &> LOG/hmm_${dataset}.log

## postprocessing:: Clean up and generate 7 col
echo "Post processing files"
perl ${perl_path}/file_ends_cleaner.pl DI/${dataset}.out DI/${dataset}.DI| perl ${perl_path}/converter_7col.pl > HMM/hmm_${dataset}.hmm

## Separate file into chromosome separated files
echo "Post processing files: Breakup files into chromosome separated files"
mkdir -p HMM/${dataset}/Chrom_Sep
awk -v Folder=HMM/${dataset} -f ../Separate_Dixon_7col_by_Chromosome.awk HMM/hmm_${dataset}.hmm

## for each chromosome file generate the final domains
#start chromosome index loop

for chrom in ${chromList[@]}; do

echo "chr${chrom}"
echo "Post processing files: Generating Domains"
perl ${perl_path}/hmm_probablity_correcter.pl HMM/${dataset}/Chrom_Sep/chr${chrom}.sep 2 0.99 ${BIN} | perl ${perl_path}/hmm-state_caller.pl ${fai_path} chr${chrom} | perl ${perl_path}/hmm-state_domains.pl > DOMAINS/${dataset}_${chrom}.domain

done

