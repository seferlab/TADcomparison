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

BIN=40000
echo "Processing size ${BIN}"

dataset=bench
mkdir -p DOMAINS/${dataset}
mkdir -p DATA/${dataset}_MrTADFinder
mkdir -p control/${dataset}_MrTADFinder

chromList=($(seq 1 100)) 

#start chromosome index loop
for chrom in ${chromList[@]}; do
python3 convert_triple_MrTADFinder.py ../benchmarks/maptest_${chrom}.txt DATA/${dataset}_MrTADFinder/maptest_${chrom}.txt
nrows=$(cat ../benchmarks/maptest_${chrom}.txt | wc -l)
python3 generate_coor.py control/${dataset}_MrTADFinder/${chrom}.file1 control/${dataset}_MrTADFinder/${chrom}.file2 ${BIN} ${nrows} chr1
/home/LiuKun/julia/julia-1.4.2/bin/julia run_MrTADFinder.jl DATA/${dataset}_MrTADFinder/maptest_${chrom}.txt control/${dataset}_MrTADFinder/${chrom}.file1 control/${dataset}_MrTADFinder/${chrom}.file2 res=1.0 1 DOMAINS/${dataset}/${dataset}_MrTADFinder.chr${chrom}
done


