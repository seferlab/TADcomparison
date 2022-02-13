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

runMrTADFinder(){
BIN=$1
echo "Processing size ${BIN}"

list=($(seq $2 $3))
chromList=($(seq 1 22)) 
#chromList[${#chromList[*]}]=X

for li in ${list[@]}; do

data=`printf "HIC%03d" $li`

mkdir -p DOMAINS/${data}
mkdir -p DATA/${data}_MrTADFinder
mkdir -p control/${data}_MrTADFinder

for chrom in ${chromList[@]}; do
python3 convert_triple_MrTADFinder.py ../Rao/${data}/${data}_50k_KR.chr${chrom} DATA/${data}_MrTADFinder/${data}_50k_KR.chr${chrom}
nrows=$(cat ../Rao/${data}/${data}_50k_KR.chr${chrom} | wc -l)
python3 generate_coor.py control/${data}_MrTADFinder/chr${chrom}.file1 control/${data}_MrTADFinder/chr${chrom}.file2 ${BIN} ${nrows} chr${chrom}
/home/LiuKun/julia/julia-1.4.2/bin/julia run_MrTADFinder.jl DATA/${data}_MrTADFinder/${data}_50k_KR.chr${chrom} control/${data}_MrTADFinder/chr${chrom}.file1 control/${data}_MrTADFinder/chr${chrom}.file2 res=1.8 ${chrom} DOMAINS/${data}/${data}_MrTADFinder.chr${chrom}
done
chrom=X
python3 convert_triple_MrTADFinder.py ../Rao/${data}/${data}_50k_KR.chr${chrom} DATA/${data}_MrTADFinder/${data}_50k_KR.chr${chrom}
nrows=$(cat ../Rao/${data}/${data}_50k_KR.chr${chrom} | wc -l)
python3 generate_coor.py control/${data}_MrTADFinder/chr${chrom}.file1 control/${data}_MrTADFinder/chr${chrom}.file2 ${BIN} ${nrows} chr${chrom}
/home/LiuKun/julia/julia-1.4.2/bin/julia run_MrTADFinder.jl DATA/${data}_MrTADFinder/${data}_50k_KR.chr${chrom} control/${data}_MrTADFinder/chr${chrom}.file1 control/${data}_MrTADFinder/chr${chrom}.file2 res=1.8 23 DOMAINS/${data}/${data}_MrTADFinder.chr${chrom}
done
}

runMrTADFinder 50000 1 29
runMrTADFinder 50000 50 56
runMrTADFinder 50000 69 74
