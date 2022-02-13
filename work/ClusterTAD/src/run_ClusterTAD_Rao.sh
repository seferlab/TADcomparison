#!/bin/bash 

checkMakeDirectory(){
        echo -e "checking directory: $1"
        if [ ! -e "$1" ]; then
                echo -e "\tmakedir $1"
                mkdir -p "$1"
        fi
}
checkMakeDirectory ../DOMAINS

runClusterTAD(){
BIN=$1
max_size=$2
echo "Processing size ${BIN}"

list=($(seq $3 $4))
chromList=($(seq 1 22)) 
chromList[${#chromList[*]}]=X

for li in ${list[@]}; do

data=`printf "HIC%03d" $li`

mkdir -p ../DOMAINS/${data}

for chrom in ${chromList[@]}; do

matlab -nodisplay -r  "filepath='../../Rao/${data}/';name='${data}_50k_KR.chr${chrom}';Res=${BIN};chromo='chr${chrom}';outputfolder_name='../DOMAINS/${data}';Max_TADsize=${max_size};",< ClusterTAD_main.m

done
done
}
runClusterTAD 50000 750000 1 29
runClusterTAD 50000 750000 50 56
runClusterTAD 50000 750000 69 74
