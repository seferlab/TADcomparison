#!/bin/bash

declare -i start_idx
declare -i end_idx
declare -i idx

PPATH=$(dirname $(readlink -f "$0"))
DPATH=${PPATH}/preprocess
mkdir -p "$DPATH/GM12878" 
	
start_idx=1
end_idx=29

for ((idx=$start_idx; idx<=$end_idx; idx++))
do
	number=$(expr $idx + 549)
	wget -P "$DPATH/GM12878" ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1551nnn/GSM1551`printf "%03d" $number`/suppl/GSM1551`printf "%03d" $number`_HIC`printf "%03d" $idx`_merged_nodups.txt.gz
done



