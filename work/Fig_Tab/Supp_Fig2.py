import numpy as np 
import pandas as pd 
import os

low_depths=[0.01,0.02,0.05]
depths=list(range(1,10))
callers=['Armatus','Arrowhead','CaTCH','CHAC','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','InsulationScore','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADbit','TADtree','TopDom']
Case_sensitive_callers=['Armatus','Arrowhead','CaTCH','''Constrained
HAC''','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','''Insulation
Score''','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADbit','TADtree','TopDom']

res=pd.DataFrame(columns = ['caller','depth', 'number'])

for depth in low_depths:
    for caller in callers:
        count = len(open('../all_TADs/bin/%s/50k_KR_downsample_ratio_%.2f.chr6'%(caller,depth), 'r').readlines())
        res=res.append({'caller':Case_sensitive_callers[callers.index(caller)],'depth':depth,'number':count},ignore_index=True)

for depth in depths:
    for caller in callers:
        count = len(open('../all_TADs/bin/%s/50k_KR_downsample_ratio_0.%d.chr6'%(caller,depth), 'r').readlines())
        res=res.append({'caller':Case_sensitive_callers[callers.index(caller)],'depth':float('%.1f' % (0.1*depth)),'number':count},ignore_index=True)
for caller in callers:
    count = len(open('../all_TADs/bin/%s/50k_KR_total.chr6'%(caller), 'r').readlines())
    res=res.append({'caller':Case_sensitive_callers[callers.index(caller)],'depth':1.0,'number':count},ignore_index=True)

#Extra Filter
res = res[res["depth"].isin([0.1, 0.5, 0.9, 1.0])] 

import matplotlib.pyplot as plt
import seaborn as sn
sn.set(style="whitegrid", color_codes=True)     #set( )设置主题，调色板更常用
sn.set_palette(sn.color_palette("Blues",13))

plt.close('all')
fig = plt.subplots(figsize=(14,9))
ax1 = plt.subplot(3,1,1)
ax2 = plt.subplot(3,1,2)
ax3 = plt.subplot(3,1,3)
# res=res[res['resolution']=='50k']
# sn.boxplot(x='caller',y='number',hue='depth',data=res,ax=ax1, showfliers=False)
sn.barplot(x='caller',y='number',hue='depth',data=res,ax=ax1)

ax1.set_ylabel('Number of TADs (#)')
ax1.set_xlabel('')
ax1.text(-0.08,1.05, 'A', transform=ax1.transAxes,fontdict = {'size': 22, 'color': 'black'})
ax1.legend(loc='right', bbox_to_anchor=(1.08, 0.4),fontsize=9)
for tick in ax1.get_xticklabels():
    tick.set_rotation(30)

res=pd.DataFrame(columns = ['caller','depth', 'size'])

for depth in low_depths:
    for caller in callers:
            chr_size=[]
            TADs=pd.read_csv('../all_TADs/bin/%s/50k_KR_downsample_ratio_%.2f.chr6'%(caller,depth),delimiter='\t', usecols=(0,1),header=None,dtype=int)
            for i in range(TADs.shape[0]):
                chr_size.append((TADs.iat[i,1]-TADs.iat[i,0]+1)/20)
            # if TADs.shape[0]==1:
            #     chr_size.append(0)
            # else:
            #     for i in range(TADs.shape[0]):
            #         chr_size.append((TADs.iat[i,1]-TADs.iat[i,0])/1000000)

            res=res.append({'caller':Case_sensitive_callers[callers.index(caller)],'depth':depth,'size':np.mean(chr_size)},ignore_index=True)

for depth in depths:
    for caller in callers:
            chr_size=[]
            TADs=pd.read_csv('../all_TADs/bin/%s/50k_KR_downsample_ratio_0.%d.chr6'%(caller,depth),delimiter='\t', usecols=(0,1),header=None,dtype=int)
            for i in range(TADs.shape[0]):
                chr_size.append((TADs.iat[i,1]-TADs.iat[i,0]+1)/20)
            # if TADs.shape[0]==1:
            #     chr_size.append(0)
            # else:
            #     for i in range(TADs.shape[0]):
            #         chr_size.append((TADs.iat[i,1]-TADs.iat[i,0])/1000000)

            res=res.append({'caller':Case_sensitive_callers[callers.index(caller)],'depth':float('%.1f' % (0.1*depth)),'size':np.mean(chr_size)},ignore_index=True)

for caller in callers:
    chr_size=[]
    TADs=pd.read_csv('../all_TADs/bin/%s/50k_KR_total.chr6'%(caller),delimiter='\t', usecols=(0,1),header=None,dtype=int)
    for i in range(TADs.shape[0]):
        chr_size.append((TADs.iat[i,1]-TADs.iat[i,0]+1)/20)
    res=res.append({'caller':Case_sensitive_callers[callers.index(caller)],'depth':1.0,'size':np.mean(chr_size)},ignore_index=True)

#Extra Filter
res = res[res["depth"].isin([0.1, 0.5, 0.9, 1.0])] 

#print(res)
#import sys
#sys.exit(1)

#res2=pd.DataFrame(columns = ['caller','depth', 'size'])

#for rind,row in res.iterrows():
#    if row["caller"] == "Matryoshka":
#        row.size = 

#import math
#res["size"] = res["size"].apply(math.log)

#import math
#for rind,row in res.iterrows():
#    row.size = math.log(row.size)
    
#print(res)
#import sys
#sys.exit(1)


sn.barplot(x='caller',y='size',hue='depth',data=res,ax=ax2)

ax2.set_ylabel('Average TAD size (Mb)')
ax2.set_xlabel('')
ax2.legend_.remove()
ax2.text(-0.08,1.05, 'B', transform=ax2.transAxes,fontdict = {'size': 22, 'color': 'black'})
for tick in ax2.get_xticklabels():
    tick.set_rotation(30)

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

res=pd.DataFrame(columns = ['caller','depth', 'number'])

from itertools import combinations
comb = list(combinations(callers, 2))

for depth in low_depths:
    for cur_comb in comb:
        ji = calculate_JI('../all_TADs/bin/%s/50k_KR_downsample_ratio_%.2f.chr6'%(cur_comb[0],depth),'../all_TADs/bin/%s/50k_KR_downsample_ratio_%.2f.chr6'%(cur_comb[1],depth))
        res=res.append({'caller':Case_sensitive_callers[callers.index(cur_comb[0])],'depth':depth,'number':ji},ignore_index=True)
        res=res.append({'caller':Case_sensitive_callers[callers.index(cur_comb[1])],'depth':depth,'number':ji},ignore_index=True)

for depth in depths:
    for cur_comb in comb:
        ji = calculate_JI('../all_TADs/bin/%s/50k_KR_downsample_ratio_0.%d.chr6'%(cur_comb[0],depth),'../all_TADs/bin/%s/50k_KR_downsample_ratio_0.%d.chr6'%(cur_comb[1],depth))
        res=res.append({'caller':Case_sensitive_callers[callers.index(cur_comb[0])],'depth':float('%.1f' % (0.1*depth)),'number':ji},ignore_index=True)
        res=res.append({'caller':Case_sensitive_callers[callers.index(cur_comb[1])],'depth':float('%.1f' % (0.1*depth)),'number':ji},ignore_index=True)
for cur_comb in comb:    
    ji = calculate_JI('../all_TADs/bin/%s/50k_KR_total.chr6'%(cur_comb[0]),'../all_TADs/bin/%s/50k_KR_total.chr6'%(cur_comb[1]))
    res=res.append({'caller':Case_sensitive_callers[callers.index(cur_comb[0])],'depth':1.0,'number':ji},ignore_index=True)
    res=res.append({'caller':Case_sensitive_callers[callers.index(cur_comb[1])],'depth':1.0,'number':ji},ignore_index=True)

#Extra Filter
res = res[res["depth"].isin([0.1, 0.5, 0.9, 1.0])] 
    
sn.boxplot(x='caller',y='number',hue='depth',data=res,ax=ax3, showfliers=False,linewidth=0.3)
ax3.set_ylabel('Jaccard Index')
ax3.set_xlabel('')
for tick in ax3.get_xticklabels():
    tick.set_rotation(30)
ax3.legend_.remove()
ax3.text(-0.08,1.05, 'C', transform=ax3.transAxes,fontdict = {'size': 22, 'color': 'black'})
plt.tight_layout()
# plt.show()
plt.savefig('Supp_Figure2.jpg',dpi=300,bbox_inches = 'tight')

