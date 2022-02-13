#!/bin/bash
PPATH=$(dirname $(readlink -f "$0"))
DPATH=${PPATH}/preprocess/GM12878

#merge raw data 
find "$DPATH" -name "*_merged_nodups.txt.gz"|xargs zcat | sort -k3,3d -k7,7d > "$DPATH/total_merged_nodups.txt"

num=$(cat $DPATH/total_merged_nodups.txt |wc -l)
num_to_split=`expr $(($num/10))`
num_to_left=`expr $(($num%10))`

#split
split -l $num_to_split $DPATH/total_merged_nodups.txt -d -a 2 $DPATH/total


#downsample
ratios=($(seq 1 9))
split_file_name=($(seq 0 9))

#for ratio in ${ratios[@]}; do
#num_downsample=`expr $(($num*$ratio/10))`
#shuf -n $num_downsample $DPATH/total_merged_nodups.txt | sort -k3,3d -k7,7d  > $DPATH/total_merged_nodups_downsample_ratio_${ratio}.txt
#done

for ratio in ${ratios[@]}; do
cat /dev/null > $DPATH/total_merged_nodups_downsample_ratio_${ratio}.tmp
num_downsample=`expr $(($num_to_split*$ratio/10))`
for split_file in ${split_file_name[@]}; do
shuf -n $num_downsample $DPATH/total0$split_file >> $DPATH/total_merged_nodups_downsample_ratio_${ratio}.tmp
done
num_downsample=`expr $(($num_to_left*$ratio/10))`
shuf -n $num_downsample $DPATH/total10 >> $DPATH/total_merged_nodups_downsample_ratio_${ratio}.tmp
cat $DPATH/total_merged_nodups_downsample_ratio_${ratio}.tmp | sort -k3,3d -k7,7d  > $DPATH/total_merged_nodups_downsample_ratio_${ratio}.txt
rm $DPATH/total_merged_nodups_downsample_ratio_${ratio}.tmp
done

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
cat /dev/null > $DPATH/total_merged_nodups_downsample_ratio_${ratio}.tmp
num_downsample=`expr $(($num_to_split*$ratio/10))`
for split_file in ${split_file_name[@]}; do
shuf -n $num_downsample $DPATH/total0$split_file >> $DPATH/total_merged_nodups_downsample_ratio_${ratio}.tmp
done
num_downsample=`expr $(($num_to_left*$ratio/10))`
shuf -n $num_downsample $DPATH/total10 >> $DPATH/total_merged_nodups_downsample_ratio_${ratio}.tmp
cat $DPATH/total_merged_nodups_downsample_ratio_${ratio}.tmp | sort -k3,3d -k7,7d  > $DPATH/total_merged_nodups_downsample_ratio_${ratio}.txt
rm $DPATH/total_merged_nodups_downsample_ratio_${ratio}.tmp
done


echo "done!"

#write your own path to juicer tools here
juicer_tool="./juicer_tools.jar"
#generate .HIC file using juicer tool, -xmx50g indicates 50g for memory which can be replaced with an appropriate value
echo "java -Xmx50g  -jar $juicer_tool pre $DPATH/total_merged_nodups.txt $DPATH/total_merged.hic hg19"
java -Xmx50g  -jar $juicer_tool pre $DPATH/total_merged_nodups.txt $DPATH/total_merged.hic hg19
ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
java -Xmx50g  -jar $juicer_tool pre $DPATH/total_merged_nodups_downsample_ratio_${ratio}.txt $DPATH/total_merged_nodups_downsample_ratio_${ratio}.hic hg19
done

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
java -Xmx50g  -jar $juicer_tool pre $DPATH/total_merged_nodups_downsample_ratio_${ratio}.txt $DPATH/total_merged_nodups_downsample_ratio_${ratio}.hic hg19
done

chrom=6
mkdir -p $DPATH/intra_KR
ratios=($(seq 1 9))
for ratio in ${ratios[@]}; do
	java -jar $juicer_tool dump observed KR $DPATH/total_merged_nodups_downsample_ratio_${ratio}.hic  $chrom $chrom BP 50000 $DPATH/intra_KR/chr${chrom}_50k_KR_downsample_ratio_${ratio}.txt -d
	python remove_nan.py $DPATH/intra_KR/chr${chrom}_50k_KR_downsample_ratio_${ratio}.txt $DPATH/intra_KR/50k_KR_downsample_ratio_${ratio}.chr${chrom}
done

ratios=(1 2 5)
for ratio in ${ratios[@]}; do
ratio=0.0${ratio}
	java -jar $juicer_tool dump observed KR $DPATH/total_merged_nodups_downsample_ratio_${ratio}.hic  $chrom $chrom BP 50000 $DPATH/intra_KR/chr${chrom}_50k_KR_downsample_ratio_${ratio}.txt -d
	python remove_nan.py $DPATH/intra_KR/chr${chrom}_50k_KR_downsample_ratio_${ratio}.txt $DPATH/intra_KR/50k_KR_downsample_ratio_${ratio}.chr${chrom}
done

resolutions=(25000 50000 100000)
for resolution in ${resolutions[@]}; do
	display_reso=`expr $(($resolution/1000))`
	java -jar $juicer_tool dump observed KR $DPATH/total_merged.hic  $chrom $chrom BP ${resolution} $DPATH/intra_KR/chr${chrom}_${resolution}_KR_total.txt -d
	python remove_nan.py $DPATH/intra_KR/chr${chrom}_${resolution}_KR_total.txt $DPATH/intra_KR/${display_reso}k_KR_total.chr${chrom}
done



