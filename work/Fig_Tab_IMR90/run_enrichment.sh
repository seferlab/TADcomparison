#!/bin/bash

mkdir -p log
mkdir -p output
nohup python3 TAD_boundary_enrichment_analysis.py CTCF marks/wgEncodeAwgTfbsSydhImr90CtcfbIggrabUniPk.narrowPeak >log/ctcf.log&
nohup python3 TAD_boundary_enrichment_analysis.py RAD21 marks/wgEncodeAwgTfbsSydhImr90Rad21IggrabUniPk.narrowPeak >log/rad21.log&
nohup python3 TAD_boundary_enrichment_analysis.py Housekeeping_genes marks/HK.txt >log/HK.log&
nohup python3 TAD_boundary_TSS_analysis.py  >log/TSS.log&
