#!/bin/bash

mkdir -p log
mkdir -p output
nohup python3 TAD_boundary_enrichment_analysis.py CTCF marks/wgEncodeAwgTfbsSydhK562CtcfbIggrabUniPk.narrowPeak >log/ctcf.log&
nohup python3 TAD_boundary_enrichment_analysis.py RAD21 marks/wgEncodeAwgTfbsSydhK562Rad21UniPk.narrowPeak >log/rad21.log&
nohup python3 TAD_boundary_enrichment_analysis.py SMC3 marks/wgEncodeAwgTfbsSydhK562Smc3ab9263IggrabUniPk.narrowPeak >log/smc3.log&
nohup python3 TAD_boundary_enrichment_analysis.py Housekeeping_genes marks/HK.txt >log/HK.log&
nohup python3 TAD_boundary_TSS_analysis.py  >log/TSS.log&
