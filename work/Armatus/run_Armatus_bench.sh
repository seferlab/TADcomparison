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
mkdir -p DATA/${dataset}_armatus

chromList=($(seq 1 100)) 

#start chromosome index loop
for chrom in ${chromList[@]}; do
cp ../benchmarks/maptest_${chrom}.txt DATA/${dataset}_armatus
mv DATA/${dataset}_armatus/maptest_${chrom}.txt DATA/${dataset}_armatus/${dataset}_armatus.chr${chrom}

gzip DATA/${dataset}_armatus/${dataset}_armatus.chr${chrom}

armatus-2.2/src/armatus -m -r ${BIN} -c ${chrom} -i ./DATA/${dataset}_armatus/${dataset}_armatus.chr${chrom}.gz -g 0.5 -o DOMAINS/${dataset}/${dataset}_armatus.chr${chrom}

done


