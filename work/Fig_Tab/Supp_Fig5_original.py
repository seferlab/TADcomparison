import numpy as np
import matplotlib.pyplot as plt
import math
import pandas
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sn

callers=['Armatus','Arrowhead','CaTCH','CHAC','CHDF','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','InsulationScore','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADbit','TADtree','TopDom']
Case_sensitive_callers=['Armatus','Arrowhead','CaTCH','''Constrained
HAC''','CHDF','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','''Insulation
Score''','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADbit','TADtree','TopDom']
# colors=['darkslategray','darkcyan','navy','c','cyan','lightgreen','lime','green','black','lightcoral','orange','darkorchid','tomato','chocolate','dodgerblue','darkred']
rainbowcolor=iter(plt.cm.Blues(np.linspace(0,1,24)))
colors=[]
for i in range(24):
    c=next(rainbowcolor)
    colors.append(c)
t = np.loadtxt('TADadjRsquared/HIC001.TADadjRsquared')
t[np.isnan(t)]=0


t=t[:,0:30]
# fig,ax = plt.subplots(1)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for i in range(len(callers)):
    lens=len(t[i,:])
    ax.plot(list(range(lens)),[i]*lens,t[i,:],c=colors[i],label=Case_sensitive_callers[i],linewidth=1)
fig.legend(loc='right', bbox_to_anchor=(1.05, 0.45),fontsize=7)
# ax.legend(loc="upper right",prop={'size': 10})
#ax.axes.get_xaxis().set_visible(False)
ax.set_zlabel('TAD$adjR^2$', {'color': 'k', 'fontsize': 12})
ax.set_xlabel('''Genomic
Distance''', {'color': 'k', 'fontsize': 10})

# plt.yticks(color='k',size=12)
plt.yticks([])
plt.xticks((0,19,29),('0','1Mb','1.5Mb'),color='k')

elev = 10
azim = 20
ax.view_init(elev, azim)

plt.tight_layout()
plt.savefig('Figure3.jpg',dpi=300, transparent=True, bbox_inches='tight')
# plt.show()