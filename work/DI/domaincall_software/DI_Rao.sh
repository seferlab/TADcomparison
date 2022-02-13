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

runDI(){
perl_path=./perl_scripts
fai_path=../hg19.fai
BIN=$1
WINDOW=$2
echo "Processing size ${BIN}"


list=($(seq $3 $4))
#chromosome index array
chromList=($(seq 1 23)) 
#chromList[${#chromList[*]}]=X


for li in ${list[@]}; do

data=`printf "HIC%03d" $li`

mkdir -p DATA/${data}_DI
mkdir -p HMM/${data}


#start chromosome index loop
for chrom in ${chromList[@]} 
do
#convert raw observed symmetric matrix to N*(N+3) matrix
if [[ chrom -eq ${chromList[${#chromList[@]}-1]} ]]; then
	python ../convert_N_N+3.py -i ../../Rao/${data}/${data}_50k_KR.chrX -o DATA/${data}_DI/${data}_DI.chr${chrom} -c chr${chrom} -b ${BIN}
else
	python ../convert_N_N+3.py -i ../../Rao/${data}/${data}_50k_KR.chr${chrom} -o DATA/${data}_DI/${data}_DI.chr${chrom} -c chr${chrom} -b ${BIN}
fi	
## Compute Directionality Index for every chromosome
echo "Computing Directionality Index chr${chrom}"
## one .matrix file for each chromosome. contains the binned genome 
perl ${perl_path}/DI_from_matrix.pl DATA/${data}_DI/${data}_DI.chr${chrom} ${BIN} ${WINDOW} ${fai_path}> DI/${data}_chr${chrom}.DI

done


## concatenate all the files together

echo "Concatenating all files together. "
cat /dev/null > DI/${data}.DI
for chrom in ${chromList[@]}; do
cat DI/${data}_chr${chrom}.DI >> DI/${data}.DI
done


## Run HMM 
echo "Running HMM for chr${chrom}"
file1="DI/${data}.DI"
file2="DI/${data}.out"
nice matlab -nodisplay -r  "inputfile='${file1}',outputfile='${file2}'",< HMM_calls.m &> LOG/hmm_${data}.log

## postprocessing:: Clean up and generate 7 col
echo "Post processing files"
perl ${perl_path}/file_ends_cleaner.pl DI/${data}.out DI/${data}.DI| perl ${perl_path}/converter_7col.pl > HMM/hmm_${data}.hmm

## Separate file into chromosome separated files
echo "Post processing files: Breakup files into chromosome separated files"
mkdir -p HMM/${data}/Chrom_Sep
awk -v Folder=HMM/${data} -f ../Separate_Dixon_7col_by_Chromosome.awk HMM/hmm_${data}.hmm

## for each chromosome file generate the final domains
#start chromosome index loop

for chrom in ${chromList[@]}; do

echo "chr${chrom}"
echo "Post processing files: Generating Domains"
perl ${perl_path}/hmm_probablity_correcter.pl HMM/${data}/Chrom_Sep/chr${chrom}.sep 2 0.99 ${BIN} | perl ${perl_path}/hmm-state_caller.pl ${fai_path} chr${chrom} | perl ${perl_path}/hmm-state_domains.pl > DOMAINS/${data}_chr${chrom}.domain

done

done
}
runDI 50000 2500000 1 29
runDI 50000 2500000 50 56
runDI 50000 2500000 69 74
