import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import pandas
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sn
sn.set(style="white", color_codes=True)     #set( )设置主题，调色板更常用
sn.set_palette(sn.color_palette("Blues",1))


callers=['Armatus','Arrowhead','CaTCH','CHAC','CHDF','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','InsulationScore','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADbit','TADtree','TopDom']
Case_sensitive_callers=['Armatus','Arrowhead','CaTCH','''Constrained
HAC''','CHDF','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','''Insulation
Score''','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADbit','TADtree','TopDom']

t = np.loadtxt('TADadjRsquared/HIC001.TADadjRsquared')
t[np.isnan(t)]=0

plt.close('all')
fig = plt.figure(figsize=(12, 5))
ax = fig.subplots(3,8)

t=t[:,0:30]
t=t.T 
t = pd.DataFrame(t,columns=callers)
t['distance'] = list(range(0,30))
t=pd.melt(t,id_vars=['distance'],value_vars=callers,var_name='callers')
for i in range(len(callers)):
    sn.lineplot(x='distance', y='value',data=t[t['callers']==callers[i]],ax=ax[int(i/8)][int(i%8)])
    ax[int(i/8)][int(i%8)].set_xlabel('')
    ax[int(i/8)][int(i%8)].set_ylabel('')
    ax[int(i/8)][int(i%8)].set_title(callers[i])
    ax[int(i/8)][int(i%8)].set_xticks((0,29))
    ax[int(i/8)][int(i%8)].set_xticklabels(('0','1.5Mb'))

    ax[int(i/8)][int(i%8)].set_ylim([0.5, 1])
    

    # ax[int(i/8)][int(i%8)].set_xticklabels((0,19,29),('0','1Mb','1.5Mb'))
# t=t[t['callers'].apply(lambda x: x in callers[0:8])]
# sn.lineplot(x='distance', y='value',hue="callers", data=t)
# ax1.legend(loc='right', bbox_to_anchor=(1.3, 0.50),fontsize=8)
sn.despine()
# plt.xticks((0,19,29),('0','1Mb','1.5Mb'),color='k')

plt.tight_layout()
plt.savefig('Supp_Figure5.jpg',dpi=300, transparent=True, bbox_inches='tight')
# plt.show()
