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

dataset=bench
mkdir -p DOMAINS/${dataset}
mkdir -p DATA/${dataset}_arrowhead

#chromosome index array
chromList=($(seq 1 100)) 

#start chromosome index loop
for chrom in ${chromList[@]}; do

echo "convert to short format with score :chr${chrom}"
python convert_arrowhead.py -i ../benchmarks/maptest_${chrom}.txt -o DATA/${dataset}_arrowhead/${dataset}_arrowhead.chr${chrom} -c chr1 -b ${BIN}


echo "gzip chr${chrom}"
gzip DATA/${dataset}_arrowhead/${dataset}_arrowhead.chr${chrom}

echo "using juicer pre command to convert chr${chrom} into gzipped file"
java -Xmx2g -jar juicer_tools.jar pre DATA/${dataset}_arrowhead/${dataset}_arrowhead.chr${chrom}.gz DATA/${dataset}_arrowhead/${dataset}_arrowhead.chr${chrom}.hic ${genomeID} -r ${BIN} -n 

echo "using juicer arrowhead command to generate domains for chr${chrom}"
java -Xmx2g -jar juicer_tools.jar arrowhead DATA/${dataset}_arrowhead/${dataset}_arrowhead.chr${chrom}.hic DOMAINS/${dataset}/${dataset}_arrowhead_${chrom}.txt -c chr1 -r ${BIN} -k NONE


done

