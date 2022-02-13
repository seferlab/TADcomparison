import numpy as np
import pandas as pd
import os

# callers=['armatus','arrowhead','CHAC','CHDF','ClusterTAD','deDoc','di','hicseg','IC-Finder','insulation','MSTD','OnTAD','SpectralTAD','TADbit','tadtree','topdom']
# Case_sensitive_callers=['Armatus','Arrowhead','''Constrained
# HAC''','CHDF','ClusterTAD','deDoc','DI','HiCseg','IC-Finder','''Insulation
# Score''','MSTD','OnTAD','SpectralTAD','TADbit','TADtree','TopDom']
# colors=['darkslategray','darkcyan','navy','c','cyan','lightgreen','lime','green','black','lightcoral','orange','darkorchid','tomato','chocolate','dodgerblue','darkred']

callers=['Armatus','Arrowhead','CaTCH','CHAC','CHDF','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','InsulationScore','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADtree','TopDom']
Case_sensitive_callers=['Armatus','Arrowhead','CaTCH','''Constrained
HAC''','CHDF','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','''Insulation
Score''','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADtree','TopDom']


c_len=len(callers)

res=np.zeros((c_len,100), dtype=float)


def calculate_TPR(file1,file2):
    # file1 is the gold standard
    count1 = len(open(file1, 'r').readlines())
    count2 = len(open(file2, 'r').readlines())
    if count1 > 0 and count2 >0:
        tads_1 = np.loadtxt(file1, delimiter='\t')
        tads_2 = np.loadtxt(file2, delimiter='\t')
        rep1=set(tads_1.flatten())
        rep2 =set(tads_2.flatten())
        intersectsize = len(set(rep1) & set(rep2))
        unionsize = len(set(rep2))
        return 1-intersectsize / unionsize
    else:
        return 0



for caller in callers:
    for idx in range(1,101):
        fdr = calculate_TPR('../benchmarks/input_seg_%d.txt'%(idx),'../all_TADs/bin/%s/benchmark.%d'%(caller,idx))
        res[callers.index(caller)][idx-1] = fdr

towrite=np.zeros((3,len(callers)))
towrite[0,:]=np.mean(res,axis=1)
towrite[1,:]=np.max(res,axis=1)
towrite[2,:]=np.min(res,axis=1)

np.savetxt('Supp_Table9_FDR.csv',towrite,fmt='%.3f',delimiter=',')

c_len=len(callers)

res=np.zeros((c_len,100), dtype=float)


def calculate_TPR(file1,file2):
    # file1 is the gold standard
    count1 = len(open(file1, 'r').readlines())
    count2 = len(open(file2, 'r').readlines())
    if count1 > 0 and count2 >0:
        tads_1 = np.loadtxt(file1, delimiter='\t')
        tads_2 = np.loadtxt(file2, delimiter='\t')
        rep1=set(tads_1.flatten())
        rep2 =set(tads_2.flatten())
        intersectsize = len(set(rep1) & set(rep2))
        unionsize = len(set(rep1))
        return intersectsize / unionsize
    else:
        return 0



for caller in callers:
    for idx in range(1,101):
        tpr = calculate_TPR('../benchmarks/input_seg_%d.txt'%(idx),'../all_TADs/bin/%s/benchmark.%d'%(caller,idx))
        res[callers.index(caller)][idx-1] = tpr

towrite=np.zeros((3,len(callers)))
towrite[0,:]=np.mean(res,axis=1)
towrite[1,:]=np.max(res,axis=1)
towrite[2,:]=np.min(res,axis=1)

np.savetxt('Supp_Table9_TPR.csv',towrite,fmt='%.3f',delimiter=',')

