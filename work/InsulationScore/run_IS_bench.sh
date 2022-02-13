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
#insulation square
square=2000000
#insulation delta span
span=800000

echo "Processing size ${BIN}"


dataset=bench

mkdir -p DOMAINS/${dataset}
mkdir -p DATA/${dataset}_insulation

#chromosome index array
chromList=($(seq 1 100)) 

#start chromosome index loop
for chrom in ${chromList[@]}; do
 
python convert_insulation.py -i ../benchmarks/maptest_${chrom}.txt -o DATA/${dataset}_insulation/${dataset}_insulation.chr${chrom}.matrix -c chr${chrom} -b ${BIN}
perl scripts/matrix2insulation.pl -i DATA/${dataset}_insulation/${dataset}_insulation.chr${chrom}.matrix -is ${square} -ids ${span} -im mean -bmoe 3 -nt 0.1 -v	
mv ${dataset}_insulation.chr${chrom}.* ./DOMAINS/${dataset}/


done

