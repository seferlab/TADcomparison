#!/bin/bash 

checkMakeDirectory(){
        echo -e "checking directory: $1"
        if [ ! -e "$1" ]; then
                echo -e "\tmakedir $1"
                mkdir -p "$1"
        fi
}
checkMakeDirectory ../DOMAINS

runSpectral(){
list=($(seq $1 $2))
chromList=($(seq 1 22)) 
chromList[${#chromList[*]}]=X

for li in ${list[@]}; do

data=`printf "HIC%03d" $li`

mkdir -p ../DOMAINS/${data}

for chrom in ${chromList[@]}; do
matlab -nodisplay -r  "file_in='../../Rao/${data}/${data}_50k_KR.chr${chrom}';file_out='../DOMAINS/${data}/${data}_Spectral.chr${chrom}'",<Spectral_run.m
done
done
}

runSpectral 1 29
runSpectral 50 56
runSpectral 69 74
