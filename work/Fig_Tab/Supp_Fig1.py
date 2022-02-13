import pandas as pd
import numpy as np
import os

def analyze_gap(file_in,binned_resolution,base_Mb):
    binned_chrom=np.zeros(binned_resolution,dtype=np.int)
    # print(binned_chrom)
    TADs = pd.read_csv(file_in, delimiter='\t', header=None)
    for i in range(TADs.shape[0]):
        start_bin=int(TADs.iat[i,0])
        end_bin=int(TADs.iat[i,1])
        for j in range(start_bin-1,end_bin):
            binned_chrom[j]=1
    gaps=pd.DataFrame(columns=['start','end'])
    previous = -1
    tmp_gap = []
    for i in range(binned_chrom.shape[0]):
        if binned_chrom[i] != previous:

            if len(tmp_gap) > 0 and previous == 0:
                gaps = gaps.append({'start': tmp_gap[0], 'end': tmp_gap[-1]}, ignore_index=True)
            tmp_gap = [i]
            previous = binned_chrom[i]
        else:
            tmp_gap.append(i)
    if len(tmp_gap) > 0 and binned_chrom[-1] == 0:
        gaps = gaps.append({'start': tmp_gap[0], 'end': tmp_gap[-1]}, ignore_index=True)
    # print(gaps)
    cur_gaps=[]
    for i in range(gaps.shape[0]):
        cur_gaps.append(((gaps.iat[i,1]-gaps.iat[i,0])+1)/base_Mb)
    return cur_gaps

dim_chr6=[6843, 3422, 1711]
resolutions=['25k','50k','100k']
num_resolution=[40,20,10]
callers=['Armatus','Arrowhead','CaTCH','CHAC','CHDF','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','InsulationScore','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADbit','TADtree','TopDom']
Case_sensitive_callers=['Armatus','Arrowhead','CaTCH','''Constrained
HAC''','CHDF','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','''Insulation
Score''','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADbit','TADtree','TopDom']

res=pd.DataFrame(columns = ['caller','resolution', 'size','num'])

for reso in resolutions:
    for caller in callers:
        if os.path.exists('../all_TADs/bin/%s/%s_KR_total.chr6'%(caller,reso)) ==True:
            gap=analyze_gap('../all_TADs/bin/%s/%s_KR_total.chr6'%(caller,reso),dim_chr6[resolutions.index(reso)],num_resolution[resolutions.index(reso)])
            to_write_size = np.mean(gap) if len(gap)>0 else 0
            # print(len(gap))
            to_write_num = len(gap)
            res=res.append({'caller':Case_sensitive_callers[callers.index(caller)],'resolution':reso,'size':to_write_size,'num':to_write_num},ignore_index=True)

import matplotlib.pyplot as plt
import seaborn as sn
sn.set(style="whitegrid", color_codes=True)     #set( )设置主题，调色板更常用
sn.set_palette(sn.color_palette("Blues",3))

plt.close('all')
fig= plt.subplots(figsize=(14,5))
ax1 = plt.subplot(2,1,1)
ax2 = plt.subplot(2,1,2)
# sn.boxplot(x='caller',y='size',hue='resolution',data=res,ax=ax1, showfliers=False)
sn.barplot(x='caller',y='num',hue='resolution',data=res,ax=ax1)

ax1.set_ylabel('Number of gaps (#)')
ax1.set_xlabel('')
ax1.legend(loc='right', bbox_to_anchor=(1.09, 0.7))
ax1.text(-0.1,1.05, 'A', transform=ax1.transAxes,fontdict = {'size': 22, 'color': 'black'})
for tick in ax1.get_xticklabels():
    tick.set_rotation(30)

# sn.boxplot(x='caller',y='num',hue='resolution',data=res,ax=ax2, showfliers=False)
sn.barplot(x='caller',y='size',hue='resolution',data=res,ax=ax2)

ax2.set_ylabel('Average gap size (Mb)')
ax2.set_xlabel('')
ax2.legend_.remove()
ax2.text(-0.1,1.05, 'B', transform=ax2.transAxes,fontdict = {'size': 22, 'color': 'black'})
# ax2.legend(loc='right', bbox_to_anchor=(1.13, 0.87))
for tick in ax2.get_xticklabels():
    tick.set_rotation(30)
plt.tight_layout()
# plt.show()
plt.savefig('Supp_Figure1.jpg',dpi=300,bbox_inches = 'tight')
# res.to_csv('test.txt')