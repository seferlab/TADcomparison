#!/bin/bash

mkdir -p TADadjRsquared
mkdir -p log

chromList2=($(seq 1 22)) 
chromList2[${#chromList2[*]}]=X

runTADadjRsquared(){
list=($(seq $1 $2))
for li in ${list[@]}; do
data=`printf "HIC%03d" $li`
nohup Rscript TADadjRsquared.r -i ../Rao/${data}/${data}_50k_KR.chr1 -d ${data} -c 22 -o TADadjRsquared/${data}.TADadjRsquared >log/${data}_TADadjRsquared.log 2>&1 &
done
}
#runTADadjRsquared 1 29
runTADadjRsquared 50 56
runTADadjRsquared 69 74
