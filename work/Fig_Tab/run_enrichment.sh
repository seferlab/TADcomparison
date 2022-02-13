#!/bin/bash

mkdir -p log
mkdir -p output
nohup python3 TAD_boundary_enrichment_analysis.py CTCF marks/wgEncodeAwgTfbsSydhGm12878Ctcfsc15914c20UniPk.narrowPeak >log/ctcf.log&
nohup python3 TAD_boundary_enrichment_analysis.py RAD21 marks/wgEncodeAwgTfbsSydhGm12878Rad21IggrabUniPk.narrowPeak >log/rad21.log&
nohup python3 TAD_boundary_enrichment_analysis.py SMC3 marks/wgEncodeAwgTfbsSydhGm12878Smc3ab9263IggmusUniPk.narrowPeak >log/smc3.log&
nohup python3 TAD_boundary_enrichment_analysis.py H3K4me3 marks/wgEncodeUwHistoneGm12878H3k4me3StdPkRep1.narrowPeak >log/H3k4me3.log&
nohup python3 TAD_boundary_enrichment_analysis.py H3K36me3 marks/wgEncodeUwHistoneGm12878H3k36me3StdPkRep1.narrowPeak >log/H3k36me3.log&
nohup python3 TAD_boundary_enrichment_analysis.py PolII marks/wgEncodeOpenChromChipGm12878Pol2Pk.narrowPeak >log/PolII.log&
nohup python3 TAD_boundary_enrichment_analysis.py Housekeeping_genes marks/HK.txt >log/HK.log&
#nohup python3 TAD_boundary_enrichment_analysis.py SINE marks/SINE.bed >log/SINE.log&
nohup python3 TAD_boundary_TSS_analysis.py  >log/TSS.log&
nohup python3 TAD_depletion_analysis.py >log/H3k9me3.log&
