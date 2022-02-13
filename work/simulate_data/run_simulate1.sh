#!/bin/bash

dataset="simulate1"
mkdir ${dataset}
Rscript SimulateData_TAD.R

list1=(409060 818120 1227180 1636240 2045301);
list2=(countMatrix interactions logs TADintervals);

for i in ${list1[@]}; do
mv simHiC_countMatrix.chr5*${i}*.txt ${dataset}/simHiC_countMatrix.chr5.${i}.txt
mv simHiC_interactions.chr5*${i}*.txt ${dataset}/simHiC_interactions.chr5.${i}.txt
mv simHiC_logs.chr5*${i}*.txt ${dataset}/simHiC_logs.chr5.${i}.txt
mv simHiC_PLOT_zoomIN.chr5*${i}*.pdf ${dataset}/simHiC_PLOT_zoomIN.chr5.${i}.pdf
mv simHiC_TADintervals.chr5*${i}*.txt ${dataset}/simHiC_TADintervals.chr5.${i}.txt
done

for i in ${list2[@]}; do
mv ${dataset}/simHiC_${i}.chr5.409060.txt ${dataset}/simHiC_${i}.chr5.0.04noise.txt
mv ${dataset}/simHiC_${i}.chr5.818120.txt ${dataset}/simHiC_${i}.chr5.0.08noise.txt
mv ${dataset}/simHiC_${i}.chr5.1227180.txt ${dataset}/simHiC_${i}.chr5.0.12noise.txt
mv ${dataset}/simHiC_${i}.chr5.1636240.txt ${dataset}/simHiC_${i}.chr5.0.16noise.txt
mv ${dataset}/simHiC_${i}.chr5.2045301.txt ${dataset}/simHiC_${i}.chr5.0.20noise.txt
done

mv ${dataset}/simHiC_PLOT_zoomIN.chr5.409060.pdf ${dataset}/simHiC_PLOT_zoomIN.chr5.0.04noise.pdf
mv ${dataset}/simHiC_PLOT_zoomIN.chr5.818120.pdf ${dataset}/simHiC_PLOT_zoomIN.chr5.0.08noise.pdf
mv ${dataset}/simHiC_PLOT_zoomIN.chr5.1227180.pdf ${dataset}/simHiC_PLOT_zoomIN.chr5.0.12noise.pdf
mv ${dataset}/simHiC_PLOT_zoomIN.chr5.1636240.pdf ${dataset}/simHiC_PLOT_zoomIN.chr5.0.16noise.pdf
mv ${dataset}/simHiC_PLOT_zoomIN.chr5.2045301.pdf ${dataset}/simHiC_PLOT_zoomIN.chr5.0.20noise.pdf

list3=(0.04 0.08 0.12 0.16 0.20);
for i in ${list3[@]}; do
python convert_simulate.py -i ${dataset}/simHiC_countMatrix.chr5.${i}noise.txt -o ${dataset}/sim_ob_${i}.chr5

done

