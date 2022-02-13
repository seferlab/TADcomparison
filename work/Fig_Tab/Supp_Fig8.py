import numpy as np
import pandas as pd
import os 
import matplotlib.pyplot as plt
import seaborn as sn
import matplotlib.patches as mpatches

callers=['Armatus','Arrowhead','CaTCH','CHAC','CHDF','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','InsulationScore','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADbit','TADtree','TopDom']
Case_sensitive_callers=['Armatus','Arrowhead','CaTCH','''Constrained
HAC''','CHDF','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','''Insulation
Score''','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADbit','TADtree','TopDom']

# rainbowcolor=iter(plt.cm.rainbow(np.linspace(0,1,24)))
# colors=[]
# for i in range(24):
#     c=next(rainbowcolor)
#     colors.append(c)

dataset=['simulate1','simulate2']

c_len=len(callers)
d_len=len(dataset)
res=np.zeros((c_len,d_len), dtype=int)
noiselist=['0.04','0.08','0.12','0.16','0.20']
noiselist=['0.04','0.08','0.12','0.16']
used_noiselist=['0.05','0.1','0.15','0.2']

res=np.zeros( (len(dataset),len(callers),len(noiselist)), dtype=int )

num_tads_res = pd.DataFrame(columns = ['data','caller','noise', 'number'])

for data in dataset:
    for ind,noise in enumerate(noiselist):
        total_number = len(open('../simulate_data/%s/simHiC_TADintervals.chr5.%snoise.txt'%(data,noise), 'r').readlines())
        total_number *= 3 # Modification for resolution
        #num_tads_res = num_tads_res.append({ 'data': data,'caller': 'Given', 'noise':used_noiselist[ind], 'number': total_number}, ignore_index=True)
        
for data in dataset:
    print('####'+data)
    for caller in callers:
        for ind,noise in enumerate(noiselist):
            if os.path.exists('../all_TADs/bin/%s/%s.%s.chr5'%(caller,data,noise)):
                total_number = len(open('../all_TADs/bin/%s/%s.%s.chr5'%(caller,data,noise), 'r').readlines())
                total_number *= 3 # Modification for resolution
                num_tads_res = num_tads_res.append({'data': data, 'caller': caller, 'noise': used_noiselist[ind], 'number': total_number}, ignore_index=True)

            else:
                num_tads_res = num_tads_res.append({'data': data, 'caller': caller, 'noise': used_noiselist[ind], 'number': 0}, ignore_index=True)

# np.savetxt('3.3_total_number.txt',res,fmt='%g',delimiter='\t')

sn.set(style="whitegrid", color_codes=True)    
colors=sn.color_palette("Blues",6)
colors=colors[-5:]
sn.set_palette(colors)

plt.close('all')
fig = plt.figure(figsize=(10,12))
# ax = fig.subplots(2,3)
# axes=[ax[0][0],ax[1][0],ax[0][1],ax[1][1],ax[0][2],ax[1][2]]
ax = fig.subplots(6,1)
axes=[ax[0],ax[3],ax[1],ax[4],ax[2],ax[5]]
# axes=[]
# for i in range(2):
#     for j in range(3):
#         axes.append(ax[i][j])

#print(num_tads_res)
#import sys
#sys.exit(1)

#3dnetmod, and pyschic missing, cunku zaten yoklar
nested_callers = ['Armatus','Arrowhead','CaTCH','deDoc', 'GMAP','Matryoshka','OnTAD','SpectralTAD','TADtree']

cnt=0
for data in dataset:
    local = num_tads_res[num_tads_res['data']==data]
    if data == "simulate2":
        local = local[local["caller"].isin(nested_callers)]        
        
    sn.despine(ax=axes[cnt])
    #sn.barplot(x='caller', y='number', hue='noise', data=num_tads_res[num_tads_res['data']==data], ax=axes[cnt])
    sn.barplot(x='caller', y='number', hue='noise', data=local, ax=axes[cnt])

    for tick in axes[cnt].get_xticklabels():
        tick.set_rotation(60)
    axes[cnt].set_ylabel('Number of TADs',size=14)
    axes[cnt].set_xlabel('')
    axes[cnt].tick_params(labelsize=10)
    for tick in axes[cnt].get_xticklabels():
        tick.set_rotation(90)
    axes[cnt].legend_.remove()
    cnt+=1

c_len=len(callers)
d_len=len(dataset)
res=np.zeros((c_len,d_len), dtype=int)
#noiselist=['0.04','0.08','0.12','0.16','0.20']

res=np.zeros( (len(dataset),len(callers),len(noiselist)), dtype=float )
TPR_res = pd.DataFrame(columns = ['data','caller','noise', 'number'])

def calculate_TPR(file1,file2):
    # file1 is the gold standard
    if os.path.exists(file2):
        count1 = len(open(file1, 'r').readlines())
        count2 = len(open(file2, 'r').readlines())
        if count1 > 0 and count2 >0:
            tads_1 = pd.read_csv(file1, sep='\t', skiprows=1, header=None)
            tads_2 = pd.read_csv(file2, sep='\t', header=None)
            rep1=set(tads_1[0])|set(tads_1[1])
            rep2=set(tads_2[0])|set(tads_2[1])
            intersectsize = len(set(rep1) & set(rep2))
            unionsize = len(set(rep1))
            return intersectsize / unionsize
        else:
            return 0
    else:
        return 0

TPR_res = pd.DataFrame(columns = ['data','caller','noise', 'number'])

for data in dataset:
    print('####'+data)
    for caller in callers:
        for ind,noise in enumerate(noiselist):
            ji_score = calculate_TPR('../simulate_data/%s/simHiC_TADintervals.chr5.%snoise.txt'%(data,noise),'../all_TADs/bin/%s/%s.%s.chr5'%(caller,data,noise))
            TPR_res = TPR_res.append({'data': data, 'caller': caller, 'noise': used_noiselist[ind], 'number': ji_score}, ignore_index=True)



for data in dataset:
    local = TPR_res[TPR_res['data']==data]
    if data == "simulate2":
        local = local[local["caller"].isin(nested_callers)]        
        
    sn.despine(ax=axes[cnt])
    sn.barplot(x='caller', y='number', hue='noise', data=local, ax=axes[cnt])
    #sn.barplot(x='caller', y='number', hue='noise', data=TPR_res[TPR_res['data']==data], ax=axes[cnt])

    axes[cnt].set_xlabel('')
    axes[cnt].set_ylabel('TPR',size=14)
    axes[cnt].tick_params(labelsize=10)
    for tick in axes[cnt].get_xticklabels():
        tick.set_rotation(90)
    axes[cnt].legend_.remove()
    cnt+=1


c_len=len(callers)
d_len=len(dataset)
res=np.zeros((c_len,d_len), dtype=int)
#noiselist=['0.04','0.08','0.12','0.16','0.20']

res=np.zeros( (len(dataset),len(callers),len(noiselist)), dtype=float )
FDR_res = pd.DataFrame(columns = ['data','caller','noise', 'number'])


def calculate_FDR(file1,file2):
    # file1 is the gold standard
    if os.path.exists(file2):
        count1 = len(open(file1, 'r').readlines())
        count2 = len(open(file2, 'r').readlines())
        if count1 > 0 and count2 >0:
            tads_1 = pd.read_csv(file1, sep='\t', skiprows=1, header=None)
            tads_2 = pd.read_csv(file2, sep='\t', header=None)
            rep1=set(tads_1[0])|set(tads_1[1])
            rep2=set(tads_2[0])|set(tads_2[1])
            intersectsize = len(set(rep1) & set(rep2))
            unionsize = len(set(rep2))
            return 1-intersectsize / unionsize
        else:
            return 0
    else:
        return 0


for data in dataset:
    print('####'+data)
    for caller in callers:
        for ind,noise in enumerate(noiselist):
            ji_score = calculate_FDR('../simulate_data/%s/simHiC_TADintervals.chr5.%snoise.txt'%(data,noise),'../all_TADs/bin/%s/%s.%s.chr5'%(caller,data,noise))
            FDR_res = FDR_res.append({'data': data, 'caller': caller, 'noise': used_noiselist[ind], 'number': ji_score}, ignore_index=True)


for data in dataset:

    local = FDR_res[FDR_res['data']==data]
    if data == "simulate2":
        local = local[local["caller"].isin(nested_callers)]        
        
    sn.despine(ax=axes[cnt])
    sn.barplot(x='caller', y='number', hue='noise', data=local, ax=axes[cnt])
    #sn.barplot(x='caller', y='number', hue='noise', data=FDR_res[FDR_res['data']==data], ax=axes[cnt])

    axes[cnt].set_xlabel('')
    axes[cnt].set_ylabel('FDR',size=14)
    axes[cnt].tick_params(labelsize=10)
    axes[cnt].legend_.remove()
    for tick in axes[cnt].get_xticklabels():
        tick.set_rotation(90)
    cnt+=1

#for i in range(5):
#    axes[i].set_xticks([])

axes[0].text(-0.2,1., 'A', transform=axes[0].transAxes,fontdict = {'size': 22, 'color': 'black'})
axes[1].text(-0.2,1., 'B', transform=axes[1].transAxes,fontdict = {'size': 22, 'color': 'black'})
axes[4].legend(loc='right', bbox_to_anchor=(1.15, 0.50))
plt.tight_layout()
plt.savefig('Supp_Figure8.jpg',dpi=300,bbox_inches = 'tight')
# plt.show()
