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
mkdir -p DATA/${dataset}

chromList=($(seq 1 100)) 

#start chromosome index loop
for chrom in ${chromList[@]}; do
python convert_homer.py -i ../benchmarks/maptest_${chrom}.txt -o DATA/${dataset}/maptest_${chrom}.txt -c chr1 -b 40000
/home/chengzj/mypython/bin/hicConvertFormat -m DATA/${dataset}/maptest_${chrom}.txt --inputFormat homer --outputFormat h5 -o DATA/${dataset}/maptest_${chrom}.h5
/home/chengzj/mypython/bin/hicFindTADs  -m DATA/${dataset}/maptest_${chrom}.h5 --outPrefix DOMAINS/${dataset}/${dataset}_HiCExplorer.chr${chrom} --numberOfProcessors 16 --correctForMultipleTesting None --minDepth 120000 --maxDepth 240000 --step 40000
done


