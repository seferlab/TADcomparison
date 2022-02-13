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
CaTCH_path=../CaTCH/DOMAINS
CHAC_path=../CHAC/DOMAINS
CHDF_path=../CHDF/DOMAINS
ClusterTAD_path=../ClusterTAD/DOMAINS
deDoc_path=../deDoc/DOMAINS
DI_path=../DI/domaincall_software/DOMAINS
EAST_path=../EAST/DOMAINS
GMAP_path=../GMAP/DOMAINS
HiCExplorer_path=../HiCExplorer/DOMAINS
HiCseg_path=../HiCseg/DOMAINS
IC_Finder_path=../IC-Finder/DOMAINS
IS_path=../InsulationScore/DOMAINS
Matryoshka_path=../Matryoshka/DOMAINS
MSTD_path=../MSTD/DOMAINS
MrTADFinder_path=../MrTADFinder/DOMAINS
OnTAD_path=../OnTAD/DOMAINS
Spectral_path=../Spectral/DOMAINS
SpectralTAD_path=../SpectralTAD/DOMAINS
TADBD_path=../TADBD/DOMAINS
TADbit_path=../TADbit/DOMAINS
TADtree_path=../TADtree/DOMAINS
TopDom_path=../TopDom/DOMAINS

simulate_data=(simulate1 simulate2)
noiselist=(0.04 0.08 0.12 0.16 0.20)

chromList2=($(seq 1 22)) 
chromList2[${#chromList2[*]}]=X

chromList4=($(seq 1 23)) 

benchList=($(seq 1 100))


checkMakeDirectory loci
checkMakeDirectory bin


#MrTADFinder
mkdir -p bin/MrTADFinder

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
python3 generate_bin_MrTADFinder.py -i ${MrTADFinder_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6 -o bin/MrTADFinder/50k_KR_downsample_ratio_${ratio}.chr6 -r 50000
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
python3 generate_bin_MrTADFinder.py -i ${MrTADFinder_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6 -o bin/MrTADFinder/50k_KR_downsample_ratio_0.${ratio}.chr6 -r 50000
done

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
python3 generate_bin_MrTADFinder.py -i ${MrTADFinder_path}/diff_reso/${display_reso}k_KR_total.chr6 -o bin/MrTADFinder/${display_reso}k_KR_total.chr6 -r ${resolution}
done

#Spectral
mkdir -p bin/Spectral

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
python3 generate_bin_Spectral.py -i ${Spectral_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6 -o bin/Spectral/50k_KR_downsample_ratio_${ratio}.chr6 
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
python3 generate_bin_Spectral.py -i ${Spectral_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6 -o bin/Spectral/50k_KR_downsample_ratio_0.${ratio}.chr6 
done

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
python3 generate_bin_Spectral.py -i ${Spectral_path}/diff_reso/${display_reso}k_KR_total.chr6 -o bin/Spectral/${display_reso}k_KR_total.chr6 
done




#MrTADFinder
extractMrTADFinder(){
mkdir -p loci/MrTADFinder
mkdir -p bin/MrTADFinder

for ((idx=$1; idx<=$2; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_MrTADFinder.py -i ${MrTADFinder_path}/${data}/${data}_MrTADFinder.chr${chrom} -o loci/MrTADFinder/${data}_MrTADFinder.chr${chrom} 
	python3 generate_bin_MrTADFinder.py -i ${MrTADFinder_path}/${data}/${data}_MrTADFinder.chr${chrom} -o bin/MrTADFinder/${data}_MrTADFinder.chr${chrom} -r 50000
	done
done
} 

#for simulated data
for data in ${simulate_data[@]}; do
for noise in ${noiselist[@]}; do
python3 generate_bin_MrTADFinder.py -i ${MrTADFinder_path}/${data}/${data}.${noise}noise -o bin/MrTADFinder/${data}.${noise}.chr5 -r 40000
done
done

#for benchmarks
for data in ${benchList[@]}; do
python3 generate_bin_MrTADFinder.py -i ${MrTADFinder_path}/bench/bench_MrTADFinder.chr${data} -o bin/MrTADFinder/benchmark.${data} -r 40000
done

#Spectral
extractSpectral(){
mkdir -p loci/Spectral
mkdir -p bin/Spectral

for ((idx=$1; idx<=$2; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_Spectral.py -i ${Spectral_path}/${data}/${data}_Spectral.chr${chrom} -o loci/Spectral/${data}_Spectral.chr${chrom} -c chr${chrom} -r 50000
	python3 generate_bin_Spectral.py -i ${Spectral_path}/${data}/${data}_Spectral.chr${chrom} -o bin/Spectral/${data}_Spectral.chr${chrom}
	done
done
} 

#for simulate data
for data in ${simulate_data[@]}; do
for noise in ${noiselist[@]}; do
python3 generate_bin_Spectral.py -i ${Spectral_path}/${data}/${data}.${noise}noise -o bin/Spectral/${data}.${noise}.chr5 
done
done

#for benchmarks
for data in ${benchList[@]}; do
python3 generate_bin_Spectral.py -i ${Spectral_path}/bench/bench_Spectral.chr${data} -o bin/Spectral/benchmark.${data}  
done

extractMrTADFinder 1 29
extractMrTADFinder 50 56
extractMrTADFinder 69 74
extractSpectral 1 29 
extractSpectral 50 56
extractSpectral 69 74

