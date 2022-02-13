#!/bin/bash 

checkMakeDirectory(){
        echo -e "checking directory: $1"
        if [ ! -e "$1" ]; then
                echo -e "\tmakedir $1"
                mkdir -p "$1"
        fi
}

checkMakeDirectory loci
checkMakeDirectory bin

chromList2=($(seq 1 22)) 
chromList2[${#chromList2[*]}]=X

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


#Armatus
extractArmatus(){
mkdir -p loci/Armatus
mkdir -p bin/Armatus

for ((idx=$1; idx<=$2; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_Armatus.py -i ${Armatus_path}/${data}/${data}_armatus.chr${chrom}.consensus.txt -o loci/Armatus/${data}_Armatus.chr${chrom}
	python3 generate_bin_Armatus.py -i ${Armatus_path}/${data}/${data}_armatus.chr${chrom}.consensus.txt -o bin/Armatus/${data}_Armatus.chr${chrom} -r 50000
	done
done
} 

#Arrowhead
extractArrowhead(){
mkdir -p loci/Arrowhead
mkdir -p bin/Arrowhead

for ((idx=$1; idx<=$2; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	arrowhead_file=${Arrowhead_path}/${data}/${data}_arrowhead.chr${chrom}/50000_blocks.bedpe
	if [ -f ${arrowhead_file} ]; then
	python3 generate_loci_Arrowhead.py -i ${arrowhead_file} -o loci/Arrowhead/${data}_Arrowhead.chr${chrom} 
	python3 generate_bin_Arrowhead.py -i ${arrowhead_file} -o bin/Arrowhead/${data}_Arrowhead.chr${chrom} -r 50000
	fi
	done
done
} 

#CaTCH 
extractCaTCH(){
mkdir -p loci/CaTCH
mkdir -p bin/CaTCH

for ((idx=$1; idx<=$2; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_CaTCH_EAST.py -i ${CaTCH_path}/${data}/${data}_CaTCH.chr${chrom} -o loci/CaTCH/${data}_CaTCH.chr${chrom} -c chr${chrom} -r 50000
	python3 generate_bin_CaTCH_EAST.py -i ${CaTCH_path}/${data}/${data}_CaTCH.chr${chrom} -o bin/CaTCH/${data}_CaTCH.chr${chrom} 
	done
done
}

#CHAC 
extractCHAC(){
mkdir -p loci/CHAC
mkdir -p bin/CHAC

for ((idx=$1; idx<=$2; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_CHAC_TADbit.py -i ${CHAC_path}/${data}/${data}.chr${chrom} -o loci/CHAC/${data}_CHAC.chr${chrom} -c chr${chrom} -r 50000
	cp ${CHAC_path}/${data}/${data}.chr${chrom} bin/CHAC/${data}_CHAC.chr${chrom} 
	done
done
} 


#CHDF
extractCHDF(){
mkdir -p loci/CHDF
mkdir -p bin/CHDF

for ((idx=$1; idx<=$2; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_CHDF.py -i ${CHDF_path}/${data}/${data}_CHDF.chr${chrom} -o loci/CHDF/${data}_CHDF.chr${chrom} -c chr${chrom} -r 50000
	cp ${CHDF_path}/${data}/${data}_CHDF.chr${chrom} bin/CHDF/${data}_CHDF.chr${chrom}
	done
done
} 

#ClusterTAD
extractClusterTAD(){
mkdir -p loci/ClusterTAD
mkdir -p bin/ClusterTAD

for ((idx=$1; idx<=$2; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_ClusterTAD.py -i ${ClusterTAD_path}/${data}/50kb/chr${chrom}/final_domain.txt -o loci/ClusterTAD/${data}_ClusterTAD.chr${chrom} -c chr${chrom} -r 50000
	python3 generate_bin_ClusterTAD.py -i ${ClusterTAD_path}/${data}/50kb/chr${chrom}/final_domain.txt -o bin/ClusterTAD/${data}_ClusterTAD.chr${chrom} 
	done
done
} 

#deDoc
extractdeDoc(){
mkdir -p loci/deDoc
mkdir -p bin/deDoc

for ((idx=$1; idx<=$2; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_deDoc.py -i ${deDoc_path}/${data}/${data}.chr${chrom}.deDoc -o loci/deDoc/${data}_deDoc.chr${chrom} -c chr${chrom} -r 50000
	python3 generate_bin_deDoc.py -i ${deDoc_path}/${data}/${data}.chr${chrom}.deDoc -o bin/deDoc/${data}_deDoc.chr${chrom} 
	done
done
} 


#DI
extractDI(){
mkdir -p loci/DI
mkdir -p bin/DI

for ((idx=$1; idx<=$2; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
		if [[ chrom -eq ${chromList2[${#chromList2[@]}-1]} ]]; then
			python3 change_chrName_DI.py -i ${DI_path}/${data}_chr23.domain -o loci/DI/${data}_DI.chrX -c chrX
			python3 generate_bin_DI.py -i ${DI_path}/${data}_chr23.domain -o bin/DI/${data}_DI.chrX -r 50000
		else
			cp ${DI_path}/${data}_chr${chrom}.domain loci/DI/${data}_DI.chr${chrom}
			python3 generate_bin_DI.py -i ${DI_path}/${data}_chr${chrom}.domain -o bin/DI/${data}_DI.chr${chrom} -r 50000
		fi
	done
done
} 

#EAST
extractEAST(){
mkdir -p loci/EAST
mkdir -p bin/EAST

for ((idx=$1; idx<=$2; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_CaTCH_EAST.py -i ${EAST_path}/${data}/${data}_EAST.chr${chrom} -o loci/EAST/${data}_EAST.chr${chrom} -c chr${chrom} -r 50000
	python3 generate_bin_CaTCH_EAST.py -i ${EAST_path}/${data}/${data}_EAST.chr${chrom} -o bin/EAST/${data}_EAST.chr${chrom} 
	done
done
} 


#GMAP
extractGMAP(){
mkdir -p loci/GMAP
mkdir -p bin/GMAP

for ((idx=$1; idx<=$2; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_GMAP.py -i ${GMAP_path}/${data}/${data}_GMAP.chr${chrom} -o loci/GMAP/${data}_GMAP.chr${chrom} -c chr${chrom} -r 50000
	python3 generate_bin_GMAP.py -i ${GMAP_path}/${data}/${data}_GMAP.chr${chrom} -o bin/GMAP/${data}_GMAP.chr${chrom} -r 50000
	done
done
} 

#HiCExplorer
extractHiCExplorer(){
mkdir -p loci/HiCExplorer
mkdir -p bin/HiCExplorer

for ((idx=$1; idx<=$2; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_HiCExplorer.py -i ${HiCExplorer_path}/${data}/${data}_HiCExplorer.chr${chrom}_domains.bed -o loci/HiCExplorer/${data}_HiCExplorer.chr${chrom} 
	python3 generate_bin_HiCExplorer.py -i ${HiCExplorer_path}/${data}/${data}_HiCExplorer.chr${chrom}_domains.bed -o bin/HiCExplorer/${data}_HiCExplorer.chr${chrom} -r 50000
	done
done
} 

#HiCseg
extractHiCseg(){
mkdir -p loci/HiCseg
mkdir -p bin/HiCseg

for ((idx=$1; idx<=$2; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_HiCseg.py -i ${HiCseg_path}/${data}/${data}.${chrom} -o loci/HiCseg/${data}_HiCseg.chr${chrom} -c chr${chrom} -r 50000
	python3 generate_bin_HiCseg.py -i ${HiCseg_path}/${data}/${data}.${chrom} -o bin/HiCseg/${data}_HiCseg.chr${chrom} 
	done
done
} 

#IC-Finder
extractIC_Finder(){
mkdir -p loci/IC-Finder
mkdir -p bin/IC-Finder
for ((idx=$1; idx<=$2; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_IC_Finder.py -i ${IC_Finder_path}/${data}/${data}_chr${chrom}_domains.txt -o loci/IC-Finder/${data}_IC-Finder.chr${chrom} -c chr${chrom} -r 50000
	cp ${IC_Finder_path}/${data}/${data}_chr${chrom}_domains.txt bin/IC-Finder/${data}_IC-Finder.chr${chrom} 
	done
done
} 

#Matryoshka
extractMatryoshka(){
mkdir -p loci/Matryoshka
mkdir -p bin/Matryoshka

for ((idx=$1; idx<=$2; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_Matryoshka.py -i ${Matryoshka_path}/${data}/${data}_Matryoshka.chr${chrom}.consensus.txt -o loci/Matryoshka/${data}_Matryoshka.chr${chrom} -r 50000
	python3 generate_bin_Matryoshka.py -i ${Matryoshka_path}/${data}/${data}_Matryoshka.chr${chrom}.consensus.txt -o bin/Matryoshka/${data}_Matryoshka.chr${chrom} -r 50000
	done
done
} 

#MSTD
extractMSTD(){
mkdir -p loci/MSTD
mkdir -p bin/MSTD

for ((idx=$1; idx<=$2; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_MSTD.py -i ${MSTD_path}/${data}/${data}.chr${chrom} -o loci/MSTD/${data}_MSTD.chr${chrom} -c chr${chrom} -r 50000
	python3 generate_bin_MSTD.py -i ${MSTD_path}/${data}/${data}.chr${chrom} -o bin/MSTD/${data}_MSTD.chr${chrom} 
	done
done
} 

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

#OnTAD
extractOnTAD(){
mkdir -p loci/OnTAD
mkdir -p bin/OnTAD

for ((idx=$1; idx<=$2; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_OnTAD.py -i ${OnTAD_path}/${data}/${data}.chr${chrom}.tad -o loci/OnTAD/${data}_OnTAD.chr${chrom} -c chr${chrom} -r 50000
	python3 generate_bin_OnTAD.py -i ${OnTAD_path}/${data}/${data}.chr${chrom}.tad -o bin/OnTAD/${data}_OnTAD.chr${chrom} 
	done
done
} 

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

#SpectralTAD
extractSpectralTAD(){
mkdir -p loci/SpectralTAD
mkdir -p bin/SpectralTAD

for ((idx=$1; idx<=$2; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_SpectralTAD.py -i ${SpectralTAD_path}/${data}/${data}.chr${chrom} -o loci/SpectralTAD/${data}_SpectralTAD.chr${chrom}
	python3 generate_bin_SpectralTAD.py -i ${SpectralTAD_path}/${data}/${data}.chr${chrom} -o bin/SpectralTAD/${data}_SpectralTAD.chr${chrom} -r 50000
	done
done
} 

#TADBD
extractTADBD(){
mkdir -p loci/TADBD
mkdir -p bin/TADBD

for ((idx=$1; idx<=$2; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_TADBD.py -i ${TADBD_path}/${data}/${data}_TADBD.chr${chrom} -o loci/TADBD/${data}_TADBD.chr${chrom} -c chr${chrom} -r 50000
	python3 generate_bin_TADBD.py -i ${TADBD_path}/${data}/${data}_TADBD.chr${chrom} -o bin/TADBD/${data}_TADBD.chr${chrom} 
	done
done
} 


#TADbit 
extractTADbit(){
mkdir -p loci/TADbit
mkdir -p bin/TADbit

for ((idx=$1; idx<=$2; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_CHAC_TADbit.py -i ${TADbit_path}/${data}/${data}.chr${chrom} -o loci/TADbit/${data}_TADbit.chr${chrom} -c chr${chrom} -r 50000
	cp ${TADbit_path}/${data}/${data}.chr${chrom} bin/TADbit/${data}_TADbit.chr${chrom} 
	done
done
} 

#TopDom
extractTopDom(){
mkdir -p loci/TopDom
mkdir -p bin/TopDom

for ((idx=$1; idx<=$2; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	python3 generate_loci_TopDom.py -i ${TopDom_path}/${data}/${data}.chr${chrom}.bed -o loci/TopDom/${data}_TopDom.chr${chrom} 
	python3 generate_bin_TopDom.py -i ${TopDom_path}/${data}/${data}.chr${chrom}.bed -o bin/TopDom/${data}_TopDom.chr${chrom}  -r 50000
	done
done
} 


extractIS(){
hg19_size=(4985 4864 3960 3821 3619 3422 3183 2927 2823 2711 2699 2677 2303 2146 2051 1806 1624 1561 1183 1260 963 1025 3106)

#IS
mkdir -p loci/InsulationScore
mkdir -p bin/InsulationScore

for ((idx=$1; idx<=$2; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	if [[ chrom -eq ${chromList2[${#chromList2[@]}-1]} ]]; then
		python3 generate_loci_IS.py -i ${IS_path}/${data}/${data}_IS_chr${chrom}.is2500001.*.boundaries -o loci/InsulationScore/${data}_InsulationScore.chr${chrom}  -c chr${chrom} -l ${hg19_size[${#hg19_size[@]}-1]} -r 50000
		python3 generate_bin_IS.py -i ${IS_path}/${data}/${data}_IS_chr${chrom}.is2500001.*.boundaries -o bin/InsulationScore/${data}_InsulationScore.chr${chrom} -l ${hg19_size[${#hg19_size[@]}-1]}
	else
		python3 generate_loci_IS.py -i ${IS_path}/${data}/${data}_IS_chr${chrom}.is2500001.*.boundaries -o loci/InsulationScore/${data}_InsulationScore.chr${chrom}  -c chr${chrom} -l ${hg19_size[${chrom}-1]} -r 50000
		python3 generate_bin_IS.py -i ${IS_path}/${data}/${data}_IS_chr${chrom}.is2500001.*.boundaries -o bin/InsulationScore/${data}_InsulationScore.chr${chrom} -l ${hg19_size[${chrom}-1]}
	fi
	done
done
} 

extractTADtree(){
hg19_tadtree=(1494 1458 1187 1145 1084 1025 953 877 846 812 809 802 690 643 614 541 486 467 353 377 287 306 930)

#TADtree
mkdir -p loci/TADtree
mkdir -p bin/TADtree

for ((idx=$1; idx<=$2; idx++)); do
	data=HIC`printf "%03d" $idx`
	for chrom in ${chromList2[@]}; do
	if [[ chrom -eq ${chromList2[${#chromList2[@]}-1]} ]]; then
		python3 generate_loci_TADtree.py -i ${TADtree_path}/${data}/chr${chrom}/N${hg19_tadtree[${#hg19_tadtree[@]}-1]}.txt -o loci/TADtree/${data}_TADtree.chr${chrom} -c chr${chrom} -r 50000
		python3 generate_bin_TADtree.py -i ${TADtree_path}/${data}/chr${chrom}/N${hg19_tadtree[${#hg19_tadtree[@]}-1]}.txt -o bin/TADtree/${data}_TADtree.chr${chrom}
	else
		python3 generate_loci_TADtree.py -i ${TADtree_path}/${data}/chr${chrom}/N${hg19_tadtree[${chrom}-1]}.txt -o loci/TADtree/${data}_TADtree.chr${chrom} -c chr${chrom} -r 50000
		python3 generate_bin_TADtree.py -i ${TADtree_path}/${data}/chr${chrom}/N${hg19_tadtree[${chrom}-1]}.txt -o bin/TADtree/${data}_TADtree.chr${chrom}
	fi
	done
done
} 


extractArmatus 50 56
extractArrowhead 50 56
extractCaTCH 50 56
extractCHAC 50 56
extractCHDF 50 56
extractClusterTAD 50 56
extractdeDoc 50 56
extractDI 50 56
extractEAST 50 56
extractGMAP 50 56
extractHiCExplorer 50 56
extractHiCseg 50 56
extractIC_Finder 50 56
extractMatryoshka 50 56
extractMSTD 50 56
extractMrTADFinder 50 56
extractOnTAD 50 56
extractSpectral 50 56
extractSpectralTAD 50 56
extractTADBD 50 56
extractTADbit 50 56
extractTopDom 50 56
extractIS 50 56
extractTADtree 50 56

extractArmatus 69 74
extractArrowhead 69 74
extractCaTCH 69 74
extractCHAC 69 74
extractCHDF 69 74
extractClusterTAD 69 74
extractdeDoc 69 74
extractDI 69 74
extractEAST 69 74
extractGMAP 69 74
extractHiCExplorer 69 74
extractHiCseg 69 74
extractIC_Finder 69 74
extractMatryoshka 69 74
extractMSTD 69 74
extractMrTADFinder 69 74
extractOnTAD 69 74
extractSpectral 69 74
extractSpectralTAD 69 74
extractTADBD 69 74
extractTADbit 69 74
extractTopDom 69 74
extractIS 69 74
extractTADtree 69 74

