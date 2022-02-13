python ../convert_N_N+3.py -i /homeb/LiuKun/zongshu/Rao/${data}/${data}_50k_KR.chr${chrom} -o DATA/${data}_matryoshka/${data}_matryoshka.chr${chrom} -c chr${chrom} -b ${BIN}

gzip DATA/${data}_matryoshka/${data}_matryoshka.chr${chrom}
/home/LiuKun/tmp/matryoshka/matryoshka-master/build/src/matryoshka -r ${BIN} -c ${chrom} -i ./DATA/${dataset}_armatus/${dataset}_armatus.chr${chrom}.gz -g 0.5 -o DOMAINS/${dataset}/${dataset}_armatus.chr${chrom}
