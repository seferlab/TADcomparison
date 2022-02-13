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

dataset=bench
mkdir -p DOMAINS/${dataset}
mkdir -p DATA/${dataset}_CaTCH

chromList=($(seq 1 100)) 

#start chromosome index loop
for chrom in ${chromList[@]}; do
python3 convert_triple_CaTCH.py  ../benchmarks/maptest_${chrom}.txt DATA/${dataset}_CaTCH/maptest_${chrom}.txt ${chrom}
gzip DATA/${dataset}_CaTCH/maptest_${chrom}.txt

Rscript synthetic_CaTCH_run.r -i DATA/${dataset}_CaTCH/maptest_${chrom}.txt.gz -o DOMAINS/${dataset}/${dataset}_CaTCH.chr${chrom}

done


