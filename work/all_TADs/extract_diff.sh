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

checkMakeDirectory bin

#Armatus
mkdir -p bin/Armatus

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
python3 generate_bin_Armatus.py -i ${Armatus_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6.consensus.txt -o bin/Armatus/50k_KR_downsample_ratio_${ratio}.chr6 -r 50000
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
python3 generate_bin_Armatus.py -i ${Armatus_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6.consensus.txt -o bin/Armatus/50k_KR_downsample_ratio_0.${ratio}.chr6 -r 50000
done

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
python3 generate_bin_Armatus.py -i ${Armatus_path}/diff_reso/${display_reso}k_KR_total.chr6.consensus.txt -o bin/Armatus/${display_reso}k_KR_total.chr6 -r ${resolution}
done


#Arrowhead
mkdir -p bin/Arrowhead

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
python3 generate_bin_Arrowhead.py -i ${Arrowhead_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6/50000_blocks.bedpe -o bin/Arrowhead/50k_KR_downsample_ratio_${ratio}.chr6 -r 50000
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
python3 generate_bin_Arrowhead.py -i ${Arrowhead_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6/50000_blocks.bedpe -o bin/Arrowhead/50k_KR_downsample_ratio_0.${ratio}.chr6 -r 50000
done

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
python3 generate_bin_Arrowhead.py -i ${Arrowhead_path}/diff_reso/${display_reso}k_KR_total.chr6/${resolution}_blocks.bedpe -o bin/Arrowhead/${display_reso}k_KR_total.chr6 -r ${resolution}
done

#CaTCH 
mkdir -p bin/CaTCH

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
python3 generate_bin_CaTCH_EAST.py -i ${CaTCH_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6 -o bin/CaTCH/50k_KR_downsample_ratio_${ratio}.chr6 
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
python3 generate_bin_CaTCH_EAST.py -i ${CaTCH_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6 -o bin/CaTCH/50k_KR_downsample_ratio_0.${ratio}.chr6
done

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
python3 generate_bin_CaTCH_EAST.py -i ${CaTCH_path}/diff_reso/${display_reso}k_KR_total.chr6 -o bin/CaTCH/${display_reso}k_KR_total.chr6 
done

#CHAC 
mkdir -p bin/CHAC

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
cp ${CHAC_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6 bin/CHAC/50k_KR_downsample_ratio_${ratio}.chr6
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
cp ${CHAC_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6 bin/CHAC/50k_KR_downsample_ratio_0.${ratio}.chr6
done

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
cp ${CHAC_path}/diff_reso/${display_reso}k_KR_total.chr6 bin/CHAC/${display_reso}k_KR_total.chr6
done

#CHDF
mkdir -p bin/CHDF

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
cp ${CHDF_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6 bin/CHDF/50k_KR_downsample_ratio_${ratio}.chr6
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
cp ${CHDF_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6 bin/CHDF/50k_KR_downsample_ratio_0.${ratio}.chr6
done

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
cp ${CHDF_path}/diff_reso/${display_reso}k_KR_total.chr6 bin/CHDF/${display_reso}k_KR_total.chr6
done

#ClusterTAD
mkdir -p bin/ClusterTAD

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
python3 generate_bin_ClusterTAD.py -i ${ClusterTAD_path}/downsample_${ratio}/50kb/chr6/final_domain.txt -o bin/ClusterTAD/50k_KR_downsample_ratio_${ratio}.chr6
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
python3 generate_bin_ClusterTAD.py -i ${ClusterTAD_path}/downsample_${ratio}/50kb/chr6/final_domain.txt -o bin/ClusterTAD/50k_KR_downsample_ratio_0.${ratio}.chr6
done

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
python3 generate_bin_ClusterTAD.py -i ${ClusterTAD_path}/diff_reso/${display_reso}kb/chr6/final_domain.txt -o bin/ClusterTAD/${display_reso}k_KR_total.chr6
done

#deDoc
mkdir -p bin/deDoc

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
python3 generate_bin_deDoc.py -i ${deDoc_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6.deDoc -o bin/deDoc/50k_KR_downsample_ratio_${ratio}.chr6
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
python3 generate_bin_deDoc.py -i ${deDoc_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6.deDoc -o bin/deDoc/50k_KR_downsample_ratio_0.${ratio}.chr6
done

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
python3 generate_bin_deDoc.py -i ${deDoc_path}/diff_reso/${display_reso}k_KR_total.chr6.deDoc -o bin/deDoc/${display_reso}k_KR_total.chr6
done

#DI
mkdir -p bin/DI

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
python3 generate_bin_DI.py -i ${DI_path}/downsample_${ratio}.chr6.domain -o bin/DI/50k_KR_downsample_ratio_${ratio}.chr6  -r 50000
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
python3 generate_bin_DI.py -i ${DI_path}/downsample_${ratio}.chr6.domain -o bin/DI/50k_KR_downsample_ratio_0.${ratio}.chr6  -r 50000
done

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
python3 generate_bin_DI.py -i ${DI_path}/${display_reso}k_KR_total.chr6.domain -o bin/DI/${display_reso}k_KR_total.chr6  -r ${resolution}
done

#EAST
mkdir -p bin/EAST

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
python3 generate_bin_CaTCH_EAST.py -i ${EAST_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6 -o bin/EAST/50k_KR_downsample_ratio_${ratio}.chr6
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
python3 generate_bin_CaTCH_EAST.py -i ${EAST_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6 -o bin/EAST/50k_KR_downsample_ratio_0.${ratio}.chr6
done

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
python3 generate_bin_CaTCH_EAST.py -i ${EAST_path}/diff_reso/${display_reso}k_KR_total.chr6 -o bin/EAST/${display_reso}k_KR_total.chr6
done

#GMAP
mkdir -p bin/GMAP

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
python3 generate_bin_GMAP.py -i ${GMAP_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6 -o bin/GMAP/50k_KR_downsample_ratio_${ratio}.chr6 -r 50000
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
python3 generate_bin_GMAP.py -i ${GMAP_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6 -o bin/GMAP/50k_KR_downsample_ratio_0.${ratio}.chr6 -r 50000
done

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
python3 generate_bin_GMAP.py -i ${GMAP_path}/diff_reso/${display_reso}k_KR_total.chr6 -o bin/GMAP/${display_reso}k_KR_total.chr6 -r ${resolution}
done

#HiCExplorer
mkdir -p bin/HiCExplorer

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
python3 generate_bin_HiCExplorer.py -i ${HiCExplorer_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6_domains.bed -o bin/HiCExplorer/50k_KR_downsample_ratio_${ratio}.chr6 -r 50000
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
python3 generate_bin_HiCExplorer.py -i ${HiCExplorer_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6_domains.bed -o bin/HiCExplorer/50k_KR_downsample_ratio_0.${ratio}.chr6 -r 50000
done

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
python3 generate_bin_HiCExplorer.py -i ${HiCExplorer_path}/diff_reso/${display_reso}k_KR_total.chr6_domains.bed -o bin/HiCExplorer/${display_reso}k_KR_total.chr6 -r ${resolution}
done

#HiCseg
mkdir -p bin/HiCseg

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
python3 generate_bin_HiCseg.py -i ${HiCseg_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6 -o bin/HiCseg/50k_KR_downsample_ratio_${ratio}.chr6 
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
python3 generate_bin_HiCseg.py -i ${HiCseg_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6 -o bin/HiCseg/50k_KR_downsample_ratio_0.${ratio}.chr6 
done

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
python3 generate_bin_HiCseg.py -i ${HiCseg_path}/diff_reso/${display_reso}k_KR_total.chr6 -o bin/HiCseg/${display_reso}k_KR_total.chr6 
done

#IC-Finder
mkdir -p bin/IC-Finder

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
cp ${IC_Finder_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6_domains.txt bin/IC-Finder/50k_KR_downsample_ratio_${ratio}.chr6 
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
cp ${IC_Finder_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6_domains.txt bin/IC-Finder/50k_KR_downsample_ratio_0.${ratio}.chr6 
done

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
cp ${IC_Finder_path}/diff_reso/${display_reso}k_KR_total.chr6_domains.txt bin/IC-Finder/${display_reso}k_KR_total.chr6 
done

#Matryoshka
mkdir -p bin/Matryoshka

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
python3 generate_bin_Matryoshka.py -i ${Matryoshka_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6.consensus.txt -o bin/Matryoshka/50k_KR_downsample_ratio_${ratio}.chr6 -r 50000
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
python3 generate_bin_Matryoshka.py -i ${Matryoshka_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6.consensus.txt -o bin/Matryoshka/50k_KR_downsample_ratio_0.${ratio}.chr6 -r 50000
done

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
python3 generate_bin_Matryoshka.py -i ${Matryoshka_path}/diff_reso/${display_reso}k_KR_total.chr6.consensus.txt -o bin/Matryoshka/${display_reso}k_KR_total.chr6 -r ${resolution} 
done

#MSTD
mkdir -p bin/MSTD

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
python3 generate_bin_MSTD.py -i ${MSTD_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6 -o bin/MSTD/50k_KR_downsample_ratio_${ratio}.chr6 
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
python3 generate_bin_MSTD.py -i ${MSTD_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6 -o bin/MSTD/50k_KR_downsample_ratio_0.${ratio}.chr6 
done

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
python3 generate_bin_MSTD.py -i ${MSTD_path}/diff_reso/${display_reso}k_KR_total.chr6 -o bin/MSTD/${display_reso}k_KR_total.chr6 
done

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


#OnTAD
mkdir -p bin/OnTAD

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
python3 generate_bin_OnTAD.py -i ${OnTAD_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6.tad -o bin/OnTAD/50k_KR_downsample_ratio_${ratio}.chr6 
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
python3 generate_bin_OnTAD.py -i ${OnTAD_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6.tad -o bin/OnTAD/50k_KR_downsample_ratio_0.${ratio}.chr6 
done

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
python3 generate_bin_OnTAD.py -i ${OnTAD_path}/diff_reso/${display_reso}k_KR_total.chr6.tad -o bin/OnTAD/${display_reso}k_KR_total.chr6 
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

#SpectralTAD
mkdir -p bin/SpectralTAD

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
python3 generate_bin_SpectralTAD.py -i ${SpectralTAD_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6 -o bin/SpectralTAD/50k_KR_downsample_ratio_${ratio}.chr6  -r 50000
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
python3 generate_bin_SpectralTAD.py -i ${SpectralTAD_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6 -o bin/SpectralTAD/50k_KR_downsample_ratio_0.${ratio}.chr6  -r 50000
done

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
python3 generate_bin_SpectralTAD.py -i ${SpectralTAD_path}/diff_reso/${display_reso}k_KR_total.chr6 -o bin/SpectralTAD/${display_reso}k_KR_total.chr6  -r ${resolution}
done

#TADBD
mkdir -p bin/TADBD

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
python3 generate_bin_TADBD.py -i ${TADBD_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6 -o bin/TADBD/50k_KR_downsample_ratio_${ratio}.chr6 
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
python3 generate_bin_TADBD.py -i ${TADBD_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6 -o bin/TADBD/50k_KR_downsample_ratio_0.${ratio}.chr6 
done

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
python3 generate_bin_TADBD.py -i ${TADBD_path}/diff_reso/${display_reso}k_KR_total.chr6 -o bin/TADBD/${display_reso}k_KR_total.chr6 
done

#TADbit 
mkdir -p bin/TADbit

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
cp ${TADbit_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6 bin/TADbit/50k_KR_downsample_ratio_${ratio}.chr6 
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
cp ${TADbit_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6 bin/TADbit/50k_KR_downsample_ratio_0.${ratio}.chr6 
done

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
cp ${TADbit_path}/diff_reso/${display_reso}k_KR_total.chr6 bin/TADbit/${display_reso}k_KR_total.chr6 
done


#TopDom
mkdir -p bin/TopDom

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
python3 generate_bin_TopDom.py -i ${TopDom_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6.bed -o bin/TopDom/50k_KR_downsample_ratio_${ratio}.chr6 -r 50000
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
python3 generate_bin_TopDom.py -i ${TopDom_path}/downsample/50k_KR_downsample_ratio_${ratio}.chr6.bed -o bin/TopDom/50k_KR_downsample_ratio_0.${ratio}.chr6 -r 50000
done

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
python3 generate_bin_TopDom.py -i ${TopDom_path}/diff_reso/${display_reso}k_KR_total.chr6.bed -o bin/TopDom/${display_reso}k_KR_total.chr6 -r ${resolution}
done


#IS
mkdir -p bin/InsulationScore

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
python3 generate_bin_IS.py -i ${IS_path}/downsample/downsample_ratio_${ratio}_chr6.is2500001.*.boundaries -o bin/InsulationScore/50k_KR_downsample_ratio_${ratio}.chr6 -l 3422
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
python3 generate_bin_IS.py -i ${IS_path}/downsample/downsample_ratio_${ratio}_chr6.is2500001.*.boundaries -o bin/InsulationScore/50k_KR_downsample_ratio_0.${ratio}.chr6 -l 3422
done

resolutions=(25000 50000 100000)
IS_chr6_len=(6843 3422 1711)
para2=0
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
para1=`expr $(($resolution*50+1))`
python3 generate_bin_IS.py -i ${IS_path}/diff_reso/${display_reso}k_chr6.is${para1}.*.boundaries -o bin/InsulationScore/${display_reso}k_KR_total.chr6 -l ${IS_chr6_len[${para2}]}
para2=`expr $(($para2+1))`
done

hg19_tadtree=(1494 1458 1187 1145 1084 1025 953 877 846 812 809 802 690 643 614 541 486 467 353 377 287 306 930)

#TADtree
mkdir -p bin/TADtree

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
python3 generate_bin_TADtree.py -i ${TADtree_path}/downsample_${ratio}/chr6/N1025.txt -o bin/TADtree/50k_KR_downsample_ratio_${ratio}.chr6 
done

ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
python3 generate_bin_TADtree.py -i ${TADtree_path}/downsample_${ratio}/chr6/N1025.txt -o bin/TADtree/50k_KR_downsample_ratio_0.${ratio}.chr6 
done

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
display_reso=`expr $(($resolution/1000))`
python3 generate_bin_TADtree.py -i ${TADtree_path}/Gm12878_${display_reso}k/chr6/N1025.txt -o bin/TADtree/${display_reso}k_KR_total.chr6 
done


