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

runAramtus(){
BIN=$1
echo "Processing size ${BIN}"

list=($(seq $2 $3))
chromList=($(seq 1 22)) 
chromList[${#chromList[*]}]=X

for li in ${list[@]}; do

data=`printf "HIC%03d" $li`

mkdir -p DOMAINS/${data}
mkdir -p DATA/${data}_armatus

for chrom in ${chromList[@]}; do
cp ../Rao/${data}/${data}_50k_KR.chr${chrom}  DATA/${data}_armatus
mv DATA/${data}_armatus/${data}_50k_KR.chr${chrom} DATA/${data}_armatus/${data}_armatus.chr${chrom}

gzip DATA/${data}_armatus/${data}_armatus.chr${chrom}

armatus-2.2/src/armatus -m -r ${BIN} -c ${chrom} -i DATA/${data}_armatus/${data}_armatus.chr${chrom}.gz -g 0.5 -o DOMAINS/${data}/${data}_armatus.chr${chrom}

done
done
}
runAramtus 50000 1 29
runAramtus 50000 50 56
runAramtus 50000 69 74
