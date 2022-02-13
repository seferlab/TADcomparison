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
BIN=40000
echo "Processing size ${BIN}"

dataset=bench
mkdir -p DOMAINS/${dataset}
mkdir -p DATA/${dataset}_Matryoshka

chromList=($(seq 1 100)) 

#start chromosome index loop
for chrom in ${chromList[@]}; do
python3 convert_N_N+3.py -i ../benchmarks/maptest_${chrom}.txt -o DATA/${dataset}_Matryoshka/maptest_${chrom}.txt  -c chr${chrom} -b ${BIN}
gzip DATA/${dataset}_Matryoshka/maptest_${chrom}.txt

/home/LiuKun/tmp/matryoshka/matryoshka-master/build/src/matryoshka -r ${BIN} -c ${chrom} -i DATA/${dataset}_Matryoshka/maptest_${chrom}.txt.gz -g 0.5 -o DOMAINS/${dataset}/${dataset}_Matryoshka.chr${chrom}

done


