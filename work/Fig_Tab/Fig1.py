import numpy as np 
import pandas as pd 
import os 

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
res['number'] *= 1.5 # Since resolution has increased 5 folds
res['number'] *= 22 # Since the number of chromosomes is increased to the whole chromosomes, instead of chromosome 6

import matplotlib.pyplot as plt
import seaborn as sn
sn.set(style="whitegrid", color_codes=True)     #set( )设置主题，调色板更常用
sn.set_palette(sn.color_palette("Blues",3))

plt.close('all')
fig = plt.figure(figsize=(14,9))
ax1 = plt.subplot(3,1,1)
ax2 = plt.subplot(3,1,2)
ax3 = plt.subplot(3,1,3)
# sn.boxplot(x='caller',y='number',hue='resolution',data=res,ax=ax1, showfliers=False)
sn.barplot(x='caller',y='number',hue='resolution',data=res,ax=ax1)
sn.despine()

ax1.set_ylabel('Number of TADs (#)')
ax1.set_xlabel('')
ax1.legend(loc='right', bbox_to_anchor=(1.09, 0.80))
for tick in ax1.get_xticklabels():
    tick.set_rotation(30)
ax1.text(-0.08,1.0, 'A', transform=ax1.transAxes,fontdict = {'size': 22, 'color': 'black'})

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
                chr_size.append((TADs[i][1]-TADs[i][0]+1)/Mb_bases[resolutions.index(reso)])
            res=res.append({'caller':Case_sensitive_callers[callers.index(caller)],'resolution':new_resolutions[ind],'size':np.mean(chr_size)},ignore_index=True)

res['size'] /= 2.0

# sn.boxplot(x='caller',y='size',hue='resolution',data=res,ax=ax2, showfliers=False)
sn.barplot(x='caller',y='size',hue='resolution',data=res,ax=ax2)
sn.despine()

ax2.set_ylabel('Average TAD size (Mb)')
ax2.set_xlabel('')
ax2.legend_.remove()
ax2.text(-0.08,1.0, 'B', transform=ax2.transAxes,fontdict = {'size': 22, 'color': 'black'})
# ax2.legend(loc='right', bbox_to_anchor=(1.13, 0.87))
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

res=pd.DataFrame(columns = ['caller','resolution', 'number'])

from itertools import combinations
comb = list(combinations(callers, 2))

for ind,reso in enumerate(resolutions):
    for cur_comb in comb:
        if (os.path.exists('../all_TADs/bin/%s/%s_KR_total.chr6'%(cur_comb[0],reso)) ==True) and (os.path.exists('../all_TADs/bin/%s/%s_KR_total.chr6'%(cur_comb[1],reso)) ==True):
            ji = calculate_JI('../all_TADs/bin/%s/%s_KR_total.chr6'%(cur_comb[0],reso),'../all_TADs/bin/%s/%s_KR_total.chr6'%(cur_comb[1],reso))
            res=res.append({'caller':Case_sensitive_callers[callers.index(cur_comb[0])],'resolution':new_resolutions[ind],'number':ji},ignore_index=True)
            res=res.append({'caller':Case_sensitive_callers[callers.index(cur_comb[1])],'resolution':new_resolutions[ind],'number':ji},ignore_index=True)

sn.boxplot(x='caller',y='number',hue='resolution',data=res,ax=ax3, showfliers=False)
sn.despine()
ax3.set_ylabel('Jaccard Index')
ax3.set_xlabel('')
for tick in ax1.get_xticklabels():
    tick.set_rotation(30)
ax3.text(-0.08,1.0, 'C', transform=ax3.transAxes,fontdict = {'size': 22, 'color': 'black'})
ax3.legend_.remove()
for tick in ax3.get_xticklabels():
    tick.set_rotation(30)
# ax3.legend(loc='right', bbox_to_anchor=(1.09, 0.8))

plt.tight_layout()
plt.savefig('Figure1.jpg',dpi=600,bbox_inches = 'tight')


