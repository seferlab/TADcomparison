import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd

def compare_TADadjR2(start_idx,end_idx,output_name):
    Case_sensitive_callers=['Armatus','Arrowhead','CaTCH','Constrained HAC','CHDF','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','InsulationScore','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADbit','TADtree','TopDom']

    dataset=list(range(start_idx,end_idx+1))
    dataset=['HIC%03d'%(x) for x in dataset]

    res=np.zeros((len(dataset),len(Case_sensitive_callers)))
    for data in dataset:
        t=pd.read_csv('TADadjRsquared/%s.TADadjRsquared'%(data),delimiter='\t',header=None)
        t=t.iloc[:,0:40]
        t=t.fillna(0)
        avg_TADadjR2=np.mean(t,axis=1)
        # print(len(avg_TADadjR2))
        res[dataset.index(data)]=avg_TADadjR2


    nres=np.zeros((3,len(Case_sensitive_callers)))
    nres[0,:]=np.mean(res,axis=0)
    nres[1,:]=np.max(res,axis=0)
    nres[2,:]=np.min(res,axis=0)
    nres = pd.DataFrame(nres,columns=Case_sensitive_callers)
    nres.to_csv(output_name,float_format='%.3f')
compare_TADadjR2(1,29,'Table5.csv')
compare_TADadjR2(50,56,'Supp_Table7.csv')
compare_TADadjR2(69,74,'Supp_Table8.csv')
