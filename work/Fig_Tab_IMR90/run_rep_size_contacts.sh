#!/bin/bash
callers=(Armatus Arrowhead CaTCH CHAC CHDF ClusterTAD deDoc DI EAST GMAP HiCExplorer HiCseg IC-Finder InsulationScore Matryoshka MrTADFinder MSTD OnTAD Spectral SpectralTAD TADBD TADbit TADtree TopDom)

mkdir -p log
mkdir -p rep_size_contacts


for caller in ${callers[@]}; do
nohup python3 rep_size_contacts.py ${caller} >log/${caller}_rep.log&
done

