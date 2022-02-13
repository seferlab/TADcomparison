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

runIS(){
BIN=$1
square=$2
span=$3
echo "Processing size ${BIN}"

list=($(seq $4 $5))
chromList=($(seq 1 22)) 
chromList[${#chromList[*]}]=X

for li in ${list[@]}; do

data=`printf "HIC%03d" $li`

mkdir -p DOMAINS/${data}
mkdir -p DATA/${data}_IS

for chrom in ${chromList[@]}; do

python convert_insulation.py -i ../Rao/${data}/${data}_50k_KR.chr${chrom} -o DATA/${data}_IS/${data}_IS_chr${chrom}.matrix -c chr${chrom} -b ${BIN}
perl scripts/matrix2insulation.pl -i DATA/${data}_IS/${data}_IS_chr${chrom}.matrix -is ${square} -ids ${span} -im mean -bmoe 3 -nt 0.1	
mv ${data}_IS_chr${chrom}.* DOMAINS/${data}/

done
done
}
runIS 50000 2500000 1000000 1 29
runIS 50000 2500000 1000000 50 56
runIS 50000 2500000 1000000 69 74

