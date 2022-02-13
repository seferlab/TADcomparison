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
echo "Processing size ${BIN}"
genomeID=hg19

datasets=(simulate1 simulate2)
for dataset in ${datasets[@]}; do


mkdir -p DOMAINS/${dataset}
mkdir -p DATA/${dataset}_arrowhead

#chromosome index array
#noise index array
noiseList=(0.04 0.08 0.12 0.16 0.20) 
#start noise index loop
for noise in ${noiseList[@]}; do

echo "convert to short format with score :chr${chrom}"
python convert_arrowhead.py -i ../simulate_data/${dataset}/sim_ob_${noise}.chr5 -o DATA/${dataset}_arrowhead/${dataset}_${noise}.chr5 -c 5 -b ${BIN}


echo "gzip chr${chrom}"
gzip DATA/${dataset}_arrowhead/${dataset}_${noise}.chr5

echo "using juicer pre command to convert chr${chrom} into gzipped file"
java -Xmx2g -jar juicer_tools.jar pre DATA/${dataset}_arrowhead/${dataset}_${noise}.chr5.gz DATA/${dataset}_arrowhead/${dataset}_${noise}.chr5.hic ${genomeID} -r ${BIN} -n 

echo "using juicer arrowhead command to generate domains for chr${chrom}"
java -Xmx2g -jar juicer_tools.jar arrowhead DATA/${dataset}_arrowhead/${dataset}_${noise}.chr5.hic DOMAINS/${dataset}/${dataset}.${noise}noise -r ${BIN} -k NONE


done

done

