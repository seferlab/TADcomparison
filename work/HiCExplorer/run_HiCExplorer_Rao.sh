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

runHiCExplorer(){
list=($(seq $1 $2))
BIN=$3
chromList=($(seq 1 22)) 
chromList[${#chromList[*]}]=X

for li in ${list[@]}; do

data=`printf "HIC%03d" $li`

mkdir -p DOMAINS/${data}
mkdir -p DATA/${data}_HiCExplorer

for chrom in ${chromList[@]}; do
para1=`expr $(($BIN*3))`
para2=`expr $(($BIN*6))`
python convert_homer.py -i ../Rao/${data}/${data}_50k_KR.chr${chrom} -o DATA/${data}_HiCExplorer/${data}.chr${chrom} -c chr${chrom} -b ${BIN}
/home/chengzj/mypython/bin/hicConvertFormat -m DATA/${data}_HiCExplorer/${data}.chr${chrom} --inputFormat homer --outputFormat h5 -o DATA/${data}_HiCExplorer/${data}.chr${chrom}.h5
/home/chengzj/mypython/bin/hicFindTADs  -m DATA/${data}_HiCExplorer/${data}.chr${chrom}.h5 --outPrefix DOMAINS/${data}/${data}_HiCExplorer.chr${chrom} --numberOfProcessors 30 --correctForMultipleTesting None --minDepth ${para1} --maxDepth ${para2} --step ${BIN}
done
done
}

runHiCExplorer 1 29 50000
runHiCExplorer 50 56 50000
runHiCExplorer 69 74 50000
