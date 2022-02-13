#!/bin/bash

#GSE63525

download_hic(){
start_idx=$1
end_idx=$2
for ((idx=$start_idx; idx<=$end_idx; idx++))
do
	number=$(expr $idx + 549)
	wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1551nnn/GSM1551`printf "%03d" $number`/suppl/GSM1551`printf "%03d" $number`_HIC`printf "%03d" $idx`.hic
done
}
download_hic 1 29
download_hic 50 56
download_hic 69 74

