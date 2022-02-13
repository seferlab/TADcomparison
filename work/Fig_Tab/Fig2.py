import numpy as np 
import pandas as pd 
import os


def calculate_JI(file1,file2):
    size1 = os.path.getsize(file1)
    size2 = os.path.getsize(file2)
    if size1==0 or size2==0:
        return 0
    tads_1 = pd.read_csv(file1, sep='\t', header=None)
    tads_2 = pd.read_csv(file2, sep='\t', header=None)
    rep1=set(tads_1[0])|set(tads_1[1])
    rep2=set(tads_2[0])|set(tads_2[1])
    intersectsize = len(set(rep1) & set(rep2))
    unionsize = len(set(rep1) | set(rep2))
    return intersectsize / unionsize

# def calculate_JI(file1,file2):
#     size1 = os.path.getsize(file1)
#     size2 = os.path.getsize(file2)
#     if size1==0 or size2==0:
#         return 0
#     tads_1 = pd.read_csv(file1, sep='\t', header=None)
#     tads_2 = pd.read_csv(file2, sep='\t', header=None)
#     rep1=set(tads_1[0]+"_"+tads_1[1].map(str))|set(tads_1[0]+"_"+tads_1[2].map(str))
#     rep2 =set(tads_2[0]+"_"+tads_2[1].map(str))|set(tads_2[0]+"_"+tads_2[2].map(str))
#     intersectsize = len(set(rep1) & set(rep2))
#     unionsize = len(set(rep1) | set(rep2))
#     return intersectsize / unionsize


#low_depths=[0.01,0.02,0.05]
#depths=list(range(1,10))
#New subset of depths
low_depths = [0.01]
depths = [1, 5, 9]

callers=['Armatus','Arrowhead','CaTCH','CHAC','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','InsulationScore','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADbit','TADtree','TopDom']
Case_sensitive_callers=['Armatus','Arrowhead','CaTCH','''Constrained
HAC''','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','''Insulation
Score''','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADbit','TADtree','TopDom']

res=pd.DataFrame(columns = ['caller','depth', 'number'])
 
for depth in low_depths:
    for caller in callers:
        ji = calculate_JI('../all_TADs/bin/%s/50k_KR_downsample_ratio_%.2f.chr6'%(caller,depth),'../all_TADs/bin/%s/50k_KR_total.chr6'%(caller) )
        res=res.append({'caller':Case_sensitive_callers[callers.index(caller)],'depth':depth,'number':ji},ignore_index=True)

for depth in depths:
    for caller in callers:
        ji = calculate_JI('../all_TADs/bin/%s/50k_KR_downsample_ratio_0.%d.chr6'%(caller,depth),'../all_TADs/bin/%s/50k_KR_total.chr6'%(caller) )
        res=res.append({'caller':Case_sensitive_callers[callers.index(caller)],'depth':float('%.1f' % (0.1*depth)),'number':ji},ignore_index=True)
for caller in callers:
    res=res.append({'caller':Case_sensitive_callers[callers.index(caller)],'depth':1.0,'number':1},ignore_index=True)

import matplotlib.pyplot as plt
import seaborn as sn
sn.set(style="whitegrid", color_codes=True)     #set( )设置主题，调色板更常用
colors=sn.color_palette("Blues",20)
colors=colors[-13:]
sn.set_palette(colors)

plt.close('all')
fig = plt.subplots(figsize=(10,3))
ax1 = plt.subplot(1,1,1)
sn.barplot(x='caller',y='number',hue='depth',data=res,ax=ax1)
sn.despine()

ax1.set_ylabel('Jaccard Index')
ax1.set_xlabel('')
for tick in ax1.get_xticklabels():
    tick.set_rotation(60)
# ax1.text(-0.12,1.05, 'a', transform=ax1.transAxes,fontdict = {'size': 22, 'color': 'black'})
ax1.legend(loc='right', bbox_to_anchor=(1.11, 0.40),fontsize=9)

plt.tight_layout()
# plt.show()
plt.savefig('Figure2.jpg',dpi=300,bbox_inches = 'tight')

