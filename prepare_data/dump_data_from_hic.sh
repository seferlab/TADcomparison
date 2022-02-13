#!/bin/bash

checkMakeDirectory(){
        echo -e "checking directory: $1"
        if [ ! -e "$1" ]; then
                echo -e "\tmakedir $1"
                mkdir -p "$1"
        fi
}

dumpdata(){
	BIN=$1
	echo "Processing size ${BIN}"

	list=($(seq $2 $3))
	chromList=($(seq 1 22)) 
	chromList[${#chromList[*]}]=X

	for li in ${list[@]}; do

	j=$((1551550+li-1))
	previous_name="GSM"$j
	latter_name=`printf "HIC%03d" $li`
	dataset=${previous_name}_${latter_name}.hic
	#echo $dataset

	mkdir -p ${latter_name}

	for chrom in ${chromList[@]}; do
	##java -jar juicer_tools.jar arrowhead /home/chengzj/juicer/work/${depth}_${dataset}/aligned/inter.hic  DOMAINS/${depth}/${dataset}/${dataset}_arrowhead.chr6 -c chr6 -r ${BIN} -k KR --ignore_sparsity
	java -jar juicer_tools.jar dump observed KR ./${dataset} chr${chrom} chr${chrom} BP 50000 ${latter_name}/${latter_name}_50k_KR.chr${chrom}_tmp -d
	python remove_nan.py ${latter_name}/${latter_name}_50k_KR.chr${chrom}_tmp ${latter_name}/${latter_name}_50k_KR.chr${chrom}

	done
	done

}

dumpdata 50000 1 29
dumpdata 50000 50 56
dumpdata 50000 69 74

