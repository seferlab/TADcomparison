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

runCaTCH(){
list=($(seq $1 $2))
chromList=($(seq 1 22)) 
chromList[${#chromList[*]}]=X

for li in ${list[@]}; do

data=`printf "HIC%03d" $li`

mkdir -p DOMAINS/${data}
mkdir -p DATA/${data}_CaTCH

for chrom in ${chromList[@]}; do
python convert_triple_CaTCH.py ../Rao/${data}/${data}_50k_KR.chr${chrom}  DATA/${data}_CaTCH/${data}_50k_KR.chr${chrom} ${chrom}
gzip DATA/${data}_CaTCH/${data}_50k_KR.chr${chrom}

Rscript CaTCH_run.r -i DATA/${data}_CaTCH/${data}_50k_KR.chr${chrom}.gz -o DOMAINS/${data}/${data}_CaTCH.chr${chrom}

done
done
}
runCaTCH 1 29
runCaTCH 50 56
runCaTCH 69 74
