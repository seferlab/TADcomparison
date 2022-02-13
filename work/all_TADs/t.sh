#!/bin/bash 

checkMakeDirectory(){
        echo -e "checking directory: $1"
        if [ ! -e "$1" ]; then
                echo -e "\tmakedir $1"
                mkdir -p "$1"
        fi
}

Armatus_path=../Armatus/DOMAINS
Arrowhead_path=../Arrowhead/DOMAINS
CHAC_path=../CHAC/DOMAINS
CHDF_path=../CHDF/DOMAINS
ClusterTAD_path=../ClusterTAD/DOMAINS
deDoc_path=../deDoc/DOMAINS
DI_path=../DI/domaincall_software/DOMAINS
HiCseg_path=../HiCseg/DOMAINS
IC_Finder_path=../IC-Finder/DOMAINS
IS_path=../InsulationScore/DOMAINS
MSTD_path=../MSTD/DOMAINS
OnTAD_path=../OnTAD/DOMAINS
SpectralTAD_path=../SpectralTAD/DOMAINS
TADbit_path=../TADbit/DOMAINS
TADtree_path=../TADtree/DOMAINS
TopDom_path=../TopDom/DOMAINS

start_idx=1
end_idx=29

simulate_data=(simulate1 simulate2)
noiselist=(0.04 0.08 0.12 0.16 0.20)

chromList2=($(seq 1 22)) 
chromList2[${#chromList2[*]}]=X

chromList4=($(seq 1 23)) 

benchList=($(seq 1 100))


checkMakeDirectory loci
checkMakeDirectory bin

#Armatus
mkdir -p loci/Armatus
mkdir -p bin/Armatus

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_Armatus.py -i ${Armatus_path}/${data}/${data}_armatus.chr${chrom}.gamma.0.5.0.txt -o loci/Armatus/${data}_Armatus.chr${chrom}
	done
done

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_bin_Armatus.py -i ${Armatus_path}/${data}/${data}_armatus.chr${chrom}.gamma.0.5.0.txt -o bin/Armatus/${data}_Armatus.chr${chrom} -r 50000
	done
done

#for simulated data
for data in ${simulate_data[@]}; do
for noise in ${noiselist[@]}; do
python3 generate_bin_Armatus.py -i ${Armatus_path}/${data}/${data}.${noise}noise.gamma.0.5.0.txt -o bin/Armatus/${data}.${noise}.chr5 -r 40000

done
done

#for benchmarks
for data in ${benchList[@]}; do
python3 generate_bin_Armatus.py -i ${Armatus_path}/bench/bench_armatus.chr${data}.gamma.0.5.0.txt -o bin/Armatus/benchmark.${data} -r 40000
done



#Arrowhead
mkdir -p loci/Arrowhead
mkdir -p bin/Arrowhead

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	arrowhead_file=${Arrowhead_path}/${data}/${data}_arrowhead.chr${chrom}/50000_blocks.bedpe
	if [ -f ${arrowhead_file} ]; then
	python3 generate_loci_Arrowhead.py -i ${arrowhead_file} -o loci/Arrowhead/${data}_Arrowhead.chr${chrom} 
	fi
	done
done

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	arrowhead_file=${Arrowhead_path}/${data}/${data}_arrowhead.chr${chrom}/50000_blocks.bedpe
	if [ -f ${arrowhead_file} ]; then
	python3 generate_bin_Arrowhead.py -i ${arrowhead_file} -o bin/Arrowhead/${data}_Arrowhead.chr${chrom} -r 50000
	fi
	done
done

#for simulate data
for data in ${simulate_data[@]}; do
for noise in ${noiselist[@]}; do
arrowhead_file=${Arrowhead_path}/${data}/${data}.${noise}noise/40000_blocks.bedpe
if [ -f ${arrowhead_file} ]; then
python3 generate_bin_Arrowhead.py -i ${arrowhead_file} -o bin/Arrowhead/${data}.${noise}.chr5 -r 40000	
fi
done
done

#for benchmarks
for data in ${benchList[@]}; do
arrowhead_file=${Arrowhead_path}/bench/bench_arrowhead_${data}.txt/40000_blocks.bedpe 
if [ -f ${arrowhead_file} ]; then
python3 generate_bin_Arrowhead.py -i ${arrowhead_file} -o bin/Arrowhead/benchmark.${data} -r 40000	
fi
done

#CHAC 
mkdir -p loci/CHAC
mkdir -p bin/CHAC

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_CHAC_TADbit.py -i ${CHAC_path}/${data}/${data}.chr${chrom} -o loci/CHAC/${data}_CHAC.chr${chrom} -c chr${chrom} -r 50000
	done
done

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	cp ${CHAC_path}/${data}/${data}.chr${chrom} bin/CHAC/${data}_CHAC.chr${chrom} 
	done
done

#for simulate data
for data in ${simulate_data[@]}; do
for noise in ${noiselist[@]}; do
cp ${CHAC_path}/${data}/${data}.${noise}noise bin/CHAC/${data}.${noise}.chr5 
done
done

#for benchmarks
for data in ${benchList[@]}; do
cp ${CHAC_path}/bench/bench_CHAC.chr${data} bin/CHAC/benchmark.${data}
done


#CHDF
mkdir -p loci/CHDF
mkdir -p bin/CHDF

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_CHDF.py -i ${CHDF_path}/${data}/${data}_CHDF.chr${chrom} -o loci/CHDF/${data}_CHDF.chr${chrom} -c chr${chrom} -r 50000
	done
done

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	cp ${CHDF_path}/${data}/${data}_CHDF.chr${chrom} bin/CHDF/${data}_CHDF.chr${chrom}
	done
done

#for simulate data
for data in ${simulate_data[@]}; do
for noise in ${noiselist[@]}; do
cp ${CHDF_path}/${data}/${data}_${noise}.chr5 bin/CHDF/${data}.${noise}.chr5
done
done

#for benchmarks
for data in ${benchList[@]}; do
cp ${CHDF_path}/bench/bench.${data} bin/CHDF/benchmark.${data}
done

#ClusterTAD
mkdir -p loci/ClusterTAD
mkdir -p bin/ClusterTAD

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_ClusterTAD.py -i ${ClusterTAD_path}/${data}/50kb/chr${chrom}/final_domain.txt -o loci/ClusterTAD/${data}_ClusterTAD.chr${chrom} -c chr${chrom} -r 50000
	done
done

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_bin_ClusterTAD.py -i ${ClusterTAD_path}/${data}/50kb/chr${chrom}/final_domain.txt -o bin/ClusterTAD/${data}_ClusterTAD.chr${chrom} 
	done
done

#for simulate data
for data in ${simulate_data[@]}; do
for noise in ${noiselist[@]}; do
python3 generate_bin_ClusterTAD.py -i ${ClusterTAD_path}/${data}/40kb/${noise}/final_domain.txt -o bin/ClusterTAD/${data}.${noise}.chr5

done
done

#for benchmarks
for data in ${benchList[@]}; do
python3 generate_bin_ClusterTAD.py -i ${ClusterTAD_path}/bench/40kb/bench${data}/final_domain.txt -o bin/ClusterTAD/benchmark.${data}
done

#deDoc
mkdir -p loci/deDoc
mkdir -p bin/deDoc

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_deDoc.py -i ${deDoc_path}/${data}/${data}.chr${chrom}.deDoc -o loci/deDoc/${data}_deDoc.chr${chrom} -c chr${chrom} -r 50000
	done
done

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_bin_deDoc.py -i ${deDoc_path}/${data}/${data}.chr${chrom}.deDoc -o bin/deDoc/${data}_deDoc.chr${chrom} 
	done
done

#for simulate data
for data in ${simulate_data[@]}; do
for noise in ${noiselist[@]}; do

python3 generate_bin_deDoc.py -i ${deDoc_path}/${data}/${data}_${noise}.chr5.deDoc -o bin/deDoc/${data}.${noise}.chr5 

done
done

#for benchmarks
for data in ${benchList[@]}; do
python3 generate_bin_deDoc.py -i ${deDoc_path}/bench/bench_${data}.txt.deDoc -o bin/deDoc/benchmark.${data}
done


#DI
mkdir -p loci/DI
mkdir -p bin/DI

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
		if [[ chrom -eq ${chromList2[${#chromList2[@]}-1]} ]]; then
			python3 change_chrName_DI.py -i ${DI_path}/${data}_chr23.domain -o loci/DI/${data}_DI.chrX -c chrX
		else
			cp ${DI_path}/${data}_chr${chrom}.domain loci/DI/${data}_DI.chr${chrom}
		fi
	done
done

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
		if [[ chrom -eq ${chromList2[${#chromList2[@]}-1]} ]]; then
			python3 generate_bin_DI.py -i ${DI_path}/${data}_chr23.domain -o bin/DI/${data}_DI.chrX -r 50000
		else
			python3 generate_bin_DI.py -i ${DI_path}/${data}_chr${chrom}.domain -o bin/DI/${data}_DI.chr${chrom} -r 50000
		fi
	done
done

#for simulate data
for data in ${simulate_data[@]}; do
for noise in ${noiselist[@]}; do
python3 generate_bin_DI.py -i ${DI_path}/${data}.${noise}noise.domain -o bin/DI/${data}.${noise}.chr5 -r 40000
done
done

#for benchmarks
for data in ${benchList[@]}; do
python3 generate_bin_DI.py -i ${DI_path}/bench_${data}.domain -o bin/DI/benchmark.${data} -r 40000
done


#HiCseg
mkdir -p loci/HiCseg
mkdir -p bin/HiCseg

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_HiCseg.py -i ${HiCseg_path}/${data}/${data}.${chrom} -o loci/HiCseg/${data}_HiCseg.chr${chrom} -c chr${chrom} -r 50000
	done
done

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_bin_HiCseg.py -i ${HiCseg_path}/${data}/${data}.${chrom} -o bin/HiCseg/${data}_HiCseg.chr${chrom} 
	done
done

#for simulate data
for data in ${simulate_data[@]}; do
for noise in ${noiselist[@]}; do
python3 generate_bin_HiCseg.py -i ${HiCseg_path}/${data}/${data}.${noise}noise -o bin/HiCseg/${data}.${noise}.chr5 
	
done
done

#for benchmarks
for data in ${benchList[@]}; do
python3 generate_bin_HiCseg.py -i ${HiCseg_path}/bench/bench_hicseg.chr${data} -o bin/HiCseg/benchmark.${data}
done

#IC-Finder
mkdir -p loci/IC-Finder
mkdir -p bin/IC-Finder
for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_IC_Finder.py -i ${IC_Finder_path}/${data}/${data}_chr${chrom}_domains.txt -o loci/IC-Finder/${data}_IC-Finder.chr${chrom} -c chr${chrom} -r 50000
	done
done

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	cp ${IC_Finder_path}/${data}/${data}_chr${chrom}_domains.txt bin/IC-Finder/${data}_IC-Finder.chr${chrom} 
	done
done

#for simulate data
for data in ${simulate_data[@]}; do
for noise in ${noiselist[@]}; do
cp ${IC_Finder_path}/${data}/${data}_${noise}_chr5_domains.txt bin/IC-Finder/${data}.${noise}.chr5 
done
done

#for benchmarks
for data in ${benchList[@]}; do
cp ${IC_Finder_path}/bench/bench_${data}_domains.txt bin/IC-Finder/benchmark.${data}
done

#MSTD
mkdir -p loci/MSTD
mkdir -p bin/MSTD

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_MSTD.py -i ${MSTD_path}/${data}/${data}.chr${chrom} -o loci/MSTD/${data}_MSTD.chr${chrom} -c chr${chrom} -r 50000
	done
done

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_bin_MSTD.py -i ${MSTD_path}/${data}/${data}.chr${chrom} -o bin/MSTD/${data}_MSTD.chr${chrom} 
	done
done

#for simulate data
for data in ${simulate_data[@]}; do
for noise in ${noiselist[@]}; do
python3 generate_bin_MSTD.py -i ${MSTD_path}/${data}/${data}_${noise}.chr5 -o bin/MSTD/${data}.${noise}.chr5 
done
done

#for benchmarks
for data in ${benchList[@]}; do
python3 generate_bin_MSTD.py -i ${MSTD_path}/bench/bench_${data}.txt -o bin/MSTD/benchmark.${data}
done


#OnTAD
mkdir -p loci/OnTAD
mkdir -p bin/OnTAD

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_OnTAD.py -i ${OnTAD_path}/${data}/${data}.chr${chrom}.tad -o loci/OnTAD/${data}_OnTAD.chr${chrom} -c chr${chrom} -r 50000
	done
done

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_bin_OnTAD.py -i ${OnTAD_path}/${data}/${data}.chr${chrom}.tad -o bin/OnTAD/${data}_OnTAD.chr${chrom} 
	done
done

#for simulate data
for data in ${simulate_data[@]}; do
for noise in ${noiselist[@]}; do
python3 generate_bin_OnTAD.py -i ${OnTAD_path}/${data}/${data}_${noise}.chr5.tad -o bin/OnTAD/${data}.${noise}.chr5 
done
done

#for benchmarks
for data in ${benchList[@]}; do
python3 generate_bin_OnTAD.py -i ${OnTAD_path}/bench/bench_${data}.txt.tad -o bin/OnTAD/benchmark.${data}
done

#SpectralTAD
mkdir -p loci/SpectralTAD
mkdir -p bin/SpectralTAD

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_SpectralTAD.py -i ${SpectralTAD_path}/${data}/${data}.chr${chrom} -o loci/SpectralTAD/${data}_SpectralTAD.chr${chrom}
	done
done

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_bin_SpectralTAD.py -i ${SpectralTAD_path}/${data}/${data}.chr${chrom} -o bin/SpectralTAD/${data}_SpectralTAD.chr${chrom} -r 50000
	done
done

#for simulate data
for data in ${simulate_data[@]}; do
for noise in ${noiselist[@]}; do

python3 generate_bin_SpectralTAD.py -i ${SpectralTAD_path}/${data}/${data}_${noise}.chr5 -o bin/SpectralTAD/${data}.${noise}.chr5 -r 40000

done
done

#for benchmarks
for data in ${benchList[@]}; do
python3 generate_bin_SpectralTAD.py -i ${SpectralTAD_path}/bench/bench_${data}.txt -o bin/SpectralTAD/benchmark.${data} -r 40000
done

#TADbit 
mkdir -p loci/TADbit
mkdir -p bin/TADbit

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_CHAC_TADbit.py -i ${TADbit_path}/${data}/${data}.chr${chrom} -o loci/TADbit/${data}_TADbit.chr${chrom} -c chr${chrom} -r 50000
	done
done

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	cp ${TADbit_path}/${data}/${data}.chr${chrom} bin/TADbit/${data}_TADbit.chr${chrom} 
	done
done

#for simulate data
for data in ${simulate_data[@]}; do
for noise in ${noiselist[@]}; do
cp ${TADbit_path}/${data}/${data}.${noise}noise bin/TADbit/${data}.${noise}.chr5 
done
done


#TopDom
mkdir -p loci/TopDom
mkdir -p bin/TopDom

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_TopDom.py -i ${TopDom_path}/${data}/${data}.chr${chrom}.bed -o loci/TopDom/${data}_TopDom.chr${chrom} 
	done
done

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_bin_TopDom.py -i ${TopDom_path}/${data}/${data}.chr${chrom}.bed -o bin/TopDom/${data}_TopDom.chr${chrom}  -r 50000
	done
done


#for simulate data
for data in ${simulate_data[@]}; do
for noise in ${noiselist[@]}; do
python3 generate_bin_TopDom.py -i ${TopDom_path}/${data}/${data}.${noise}noise.bed -o bin/TopDom/${data}.${noise}.chr5 -r 40000
		
done
done

#for benchmarks
for data in ${benchList[@]}; do
python3 generate_bin_TopDom.py -i ${TopDom_path}/bench/bench_topdom.chr${data}.bed -o bin/TopDom/benchmark.${data} -r 40000
done

hg19_size=(4985 4864 3960 3821 3619 3422 3183 2927 2823 2711 2699 2677 2303 2146 2051 1806 1624 1561 1183 1260 963 1025 3106)

#IS
mkdir -p loci/InsulationScore
mkdir -p bin/InsulationScore

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	if [[ chrom -eq ${chromList2[${#chromList2[@]}-1]} ]]; then
		python3 generate_loci_IS.py -i ${IS_path}/${data}/${data}_IS_chr${chrom}.is2500001.*.boundaries -o loci/InsulationScore/${data}_InsulationScore.chr${chrom}  -c chr${chrom} -l ${hg19_size[${#hg19_size[@]}-1]} -r 50000
	else
		python3 generate_loci_IS.py -i ${IS_path}/${data}/${data}_IS_chr${chrom}.is2500001.*.boundaries -o loci/InsulationScore/${data}_InsulationScore.chr${chrom}  -c chr${chrom} -l ${hg19_size[${chrom}-1]} -r 50000
	fi
	done
done

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	if [[ chrom -eq ${chromList2[${#chromList2[@]}-1]} ]]; then
		python3 generate_bin_IS.py -i ${IS_path}/${data}/${data}_IS_chr${chrom}.is2500001.*.boundaries -o bin/InsulationScore/${data}_InsulationScore.chr${chrom} -l ${hg19_size[${#hg19_size[@]}-1]}
	else
		python3 generate_bin_IS.py -i ${IS_path}/${data}/${data}_IS_chr${chrom}.is2500001.*.boundaries -o bin/InsulationScore/${data}_InsulationScore.chr${chrom} -l ${hg19_size[${chrom}-1]}
	fi
	done
done

#for simulate data
for data in ${simulate_data[@]}; do
for noise in ${noiselist[@]}; do

count=$(cat ../simulate_data/${data}/simHiC_countMatrix.chr5.${noise}noise.txt |wc -l)
python3 generate_bin_IS.py -i ${IS_path}/${data}/${data}_${noise}.chr5.is2000001.*.boundaries -o bin/InsulationScore/${data}.${noise}.chr5 -l ${count}
		
done
done


#for benchmarks
for data in ${benchList[@]}; do
python3 generate_bin_IS.py -i ${IS_path}/bench/bench_insulation.chr${data}.is2000001.*.boundaries -o bin/InsulationScore/benchmark.${data} -l 500
done

hg19_tadtree=(1494 1458 1187 1145 1084 1025 953 877 846 812 809 802 690 643 614 541 486 467 353 377 287 306 930)

#TADtree
mkdir -p loci/TADtree
mkdir -p bin/TADtree

for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	if [[ chrom -eq ${chromList2[${#chromList2[@]}-1]} ]]; then
		python3 generate_loci_TADtree.py -i ${TADtree_path}/${data}/chr${chrom}/N${hg19_tadtree[${#hg19_tadtree[@]}-1]}.txt -o loci/TADtree/${data}_TADtree.chr${chrom} -c chr${chrom} -r 50000
	else
		python3 generate_loci_TADtree.py -i ${TADtree_path}/${data}/chr${chrom}/N${hg19_tadtree[${chrom}-1]}.txt -o loci/TADtree/${data}_TADtree.chr${chrom} -c chr${chrom} -r 50000
	fi
	
	done
done


for ((idx=$start_idx; idx<=$end_idx; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	if [[ chrom -eq ${chromList2[${#chromList2[@]}-1]} ]]; then
		python3 generate_bin_TADtree.py -i ${TADtree_path}/${data}/chr${chrom}/N${hg19_tadtree[${#hg19_tadtree[@]}-1]}.txt -o bin/TADtree/${data}_TADtree.chr${chrom}
	else
		python3 generate_bin_TADtree.py -i ${TADtree_path}/${data}/chr${chrom}/N${hg19_tadtree[${chrom}-1]}.txt -o bin/TADtree/${data}_TADtree.chr${chrom}
	fi
	
	done
done


#for simulate data
for data in ${simulate_data[@]}; do

python3 generate_bin_TADtree.py -i ${TADtree_path}/${data}/chr1/N1084.txt -o bin/TADtree/${data}.0.04.chr5
python3 generate_bin_TADtree.py -i ${TADtree_path}/${data}/chr2/N1084.txt -o bin/TADtree/${data}.0.08.chr5
python3 generate_bin_TADtree.py -i ${TADtree_path}/${data}/chr3/N1084.txt -o bin/TADtree/${data}.0.12.chr5
python3 generate_bin_TADtree.py -i ${TADtree_path}/${data}/chr4/N1084.txt -o bin/TADtree/${data}.0.16.chr5
python3 generate_bin_TADtree.py -i ${TADtree_path}/${data}/chr5/N1084.txt -o bin/TADtree/${data}.0.20.chr5

done

#for benchmarks
for data in ${benchList[@]}; do
python3 generate_bin_TADtree.py -i ${TADtree_path}/benchmarks/bench${data}/N49.txt -o bin/TADtree/benchmark.${data}
done



