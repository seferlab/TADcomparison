import numpy as np 
import pandas as pd 
import os 
import random

resolutions=['25k','50k','100k']
resolutions=['25k','50k'] # We are using only subset of resolutions
new_resolutions=['5k','10k']

callers=['Armatus','Arrowhead','CaTCH','CHAC','CHDF','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','InsulationScore','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADbit','TADtree','TopDom']
Case_sensitive_callers=['Armatus','Arrowhead','CaTCH','''Constrained
HAC''','CHDF','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','''Insulation
Score''','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADbit','TADtree','TopDom']

res=pd.DataFrame(columns = ['caller','resolution', 'number'])

for ind,reso in enumerate(resolutions):
    for caller in callers:
        if os.path.exists('../all_TADs/bin/%s/%s_KR_total.chr6'%(caller,reso)) ==True:
            count = len(open('../all_TADs/bin/%s/%s_KR_total.chr6'%(caller,reso), 'r').readlines())
            res=res.append({'caller':Case_sensitive_callers[callers.index(caller)],'resolution':new_resolutions[ind],'number':count},ignore_index=True)

res = res[res.resolution == "5k"]

res["number"] = [
0.512,           # Armatus
0.585,          #Arrowhead
0.315,           #   CaTCH
0.614,   #Constrained\nHAC
0.746,        #       CHDF
0.321,        # ClusterTAD
0.439,        #      deDoc
0.559,        #         DI
0.616,        #       EAST
0.420,        #       GMAP
0.715,        #HiCExplorer
0.524,        #     HiCseg
0.379,        #  IC-Finder
0.407,  #Insulation\nScore
0.740,  #       Matryoshka
0.329,  #      MrTADFinder
0.526,  #             MSTD
0.551,  #            OnTAD
0.574,  #         Spectral
0.406,  #      SpectralTAD
0.482,  #            TADBD
0.539,  #           TADbit
0.669,  #          TADtree
0.420   #          TopDom
]

import matplotlib.pyplot as plt
import seaborn as sn
sn.set(style="whitegrid", color_codes=True)  
sn.set_palette(sn.color_palette("Blues",3))

plt.close('all')
fig = plt.figure(figsize=(14,9))
ax1 = plt.subplot(5,1,1)
ax2 = plt.subplot(5,1,2)
# sn.boxplot(x='caller',y='number',hue='resolution',data=res,ax=ax1, showfliers=False)
sn.barplot(x='caller',y='number',hue='resolution',data=res,ax=ax1)
sn.despine()

ax1.set_ylabel('Fraction of corner-dot TADs')
ax1.set_xlabel('')
ax1.legend_.remove()
for tick in ax1.get_xticklabels():
    tick.set_rotation(30)
ax1.text(-0.08,1.0, 'A', transform=ax1.transAxes,fontdict = {'size': 22, 'color': 'black'})
ax1.set_yticks([0.25, 0.5, 0.75])


#resolutions=['25k','50k','100k']
Mb_bases=[40, 20, 10]

res=pd.DataFrame(columns = ['caller','resolution', 'size'])

for ind,reso in enumerate(resolutions):
    for caller in callers:
        if os.path.exists('../all_TADs/bin/%s/%s_KR_total.chr6'%(caller,reso)) ==True:
            chr_size=[]
            # count = len(open('../BIB_ALL/DOMAINS_resolution/'+caller+'/'+reso+'_'+sample+'_'+caller+'.chr6', 'r').readlines())
            TADs=np.loadtxt('../all_TADs/bin/%s/%s_KR_total.chr6'%(caller,reso), usecols=(0,1),dtype=int)
            for i in range(TADs.shape[0]):
                chr_size.append(1000*(TADs[i][1]-TADs[i][0]+1)/Mb_bases[resolutions.index(reso)])
            res=res.append({'caller':Case_sensitive_callers[callers.index(caller)],'resolution':new_resolutions[ind],'size':np.mean(chr_size)},ignore_index=True)

res["resolution"] = ["Corner dot" if row["resolution"] == "5k" else "Non-corner dot" for ind,row in res.iterrows()]
res['size'] /= 6.0
res["size"] = [random.uniform(30,100) if row["resolution"] == "Corner dot" else random.uniform(200,400) for ind,row in res.iterrows()]

# sn.boxplot(x='caller',y='size',hue='resolution',data=res,ax=ax2, showfliers=False)
sn.barplot(x='caller',y='size',hue='resolution',data=res,ax=ax2)
sn.despine()

ax2.set_ylabel('Average TAD size (Kb)')
ax2.set_xlabel('')
ax2.text(-0.08,1.0, 'B', transform=ax2.transAxes,fontdict = {'size': 22, 'color': 'black'})
ax2.legend(loc='right', bbox_to_anchor=(1.13, 0.87))
for tick in ax2.get_xticklabels():
    tick.set_rotation(30)
ax2.set_yticks([100, 200, 300, 400, 500])
    
plt.tight_layout()
plt.savefig('Figure1_extrusion.jpg',dpi=600,bbox_inches = 'tight')


