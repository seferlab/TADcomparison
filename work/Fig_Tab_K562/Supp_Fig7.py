import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from matplotlib.gridspec import GridSpec
from matplotlib.pyplot import MultipleLocator

plt.close('all')
# fig, ax1 = plt.subplots()
fig = plt.figure(figsize=(12,8))
gs = GridSpec(4,1)
ax2 = plt.subplot(gs[0, 0])
ax3 = plt.subplot(gs[1, 0])
ax4 = plt.subplot(gs[2, 0])
ax5 = plt.subplot(gs[3, 0])


callers=['Armatus','Arrowhead','CaTCH','CHAC','CHDF','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','InsulationScore','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADbit','TADtree','TopDom']
Case_sensitive_callers=['Armatus','Arrowhead','CaTCH','Constrained HAC','CHDF','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','InsulationScore','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADbit','TADtree','TopDom']

# colors=['darkslategray','darkcyan','navy','c','cyan','lightgreen','lime','green','black','lightcoral','orange','darkorchid','tomato','chocolate','dodgerblue','darkred']
# rainbowcolor=iter(plt.cm.rainbow(np.linspace(0,1,24)))
# colors=[]
# for i in range(24):
#     c=next(rainbowcolor)
#     colors.append(c)

chrs=list(range(1,23))
chrs = [str(x) for x in chrs]
chrs.append('X')


# dataset=list(range(1,30))
# dataset=['HIC%03d'%(x) for x in dataset]
dataset=['HIC071','HIC072','HIC073','HIC074']

tab='rep'

for caller in callers:
    TAD_dict = {}
    for data in dataset:
        for chr in chrs:
            if os.path.exists('../all_TADs/bin/%s/%s_%s.chr%s'%(caller,data,caller,chr))==True:
                current_TAD = pd.read_csv('../all_TADs/bin/%s/%s_%s.chr%s'%(caller,data,caller,chr), sep='\t', header=None)
                current_TAD = current_TAD.drop_duplicates()
                for index in range(0, current_TAD.shape[0]):
                    if ('chr' + chr,current_TAD.iat[index, 0], current_TAD.iat[index, 1]) in TAD_dict.keys():
                        TAD_dict[('chr' + chr, current_TAD.iat[index, 0], current_TAD.iat[index, 1])] += 1
                    else:
                        TAD_dict[('chr' + chr, current_TAD.iat[index, 0], current_TAD.iat[index, 1])] = 1
    TAD_dict = [[x[0], x[1], x[2], TAD_dict[x]] for x in TAD_dict.keys()]
    TAD_dict = pd.DataFrame(TAD_dict)
    TAD_dict = TAD_dict.sort_values(by=[0, 3], ascending=True)
    TAD_dict.to_csv('%s_%s.txt'%(tab,caller), sep='\t', header=None, index=None)



res_ratio = pd.DataFrame(columns=('callers', 'rep', 'value'))
for caller in callers:
    repro = pd.read_csv( '%s_%s.txt'%(tab,caller), sep='\t', header=None, usecols=[3])
    row_repro = repro.shape[0]
    unique_count = repro.nunique()[3]
    for i in range(1, unique_count + 1):
        #current_cnt = repro[3].apply(lambda x: 1 if x == i else 0).sum()
        current_cnt=repro[repro[3]==i].shape[0]
        res_ratio = res_ratio.append([{'callers':caller,'rep':i,'value':(current_cnt / row_repro)*100}], ignore_index=True)


res_ratio=res_ratio.set_index(['callers','rep'])['value'].unstack()
res_ratio=res_ratio.reindex(index=callers)
res_ratio.plot(kind='bar',stacked = True,ax=ax2,alpha =0.6)

ax2.set_ylabel('Ratio(%)', size=14)
# ax2.tick_params(labelsize=14)

# ax2.legend(loc="upper right",prop={'size': 13})
# ax2.legend(loc='right', bbox_to_anchor=(1.09, 0.8))
ax2.set_xlabel('')
ax2.legend_.remove()
for tick in ax2.get_xticklabels():
    tick.set_rotation(60)

ax2.text(-0.1,1.05, 'A', transform=ax2.transAxes,fontdict = {'size': 22, 'color': 'black'})



'''
import itertools

tmp1=['HIC003','HIC005','HIC014'] 

bio1={1:tmp1,2:['HIC020'],3:['HIC022','HIC023'],4:['HIC025','HIC026'],5:['HIC027']}
li=list(range(1,6))

tmp2=itertools.combinations(tmp1, 2)  
tmp3=np.zeros(shape=(0,2))
for i in tmp2:
    tmp3=np.row_stack((tmp3,[i[0],i[1]]))

c1_tec_rep=np.array([['HIC022','HIC023'],
                    ['HIC025','HIC026']])
c1_tec_rep=np.row_stack((tmp3, c1_tec_rep))


c_len=len(callers)

JI_c1_bio=np.zeros((31*23,c_len))
JI_c1_tec=np.zeros((c1_tec_rep.shape[0]*23,c_len))
'''
chrs1=list(range(1,23))
chrs2 = [str(x) for x in chrs1]
chrs2.append('X')

def calculate_JI(file1,file2):
    if os.path.exists(file1)==False or os.path.exists(file2)==False:
        return 0
    count1 = len(open(file1, 'r').readlines())
    count2 = len(open(file2, 'r').readlines())
    if count1 > 0 and count2 >0:
        tads_1 = pd.read_csv(file1, sep='\t', header=None)
        tads_2 = pd.read_csv(file2, sep='\t', header=None)
        rep1=set(tads_1[0])|set(tads_1[1])
        rep2=set(tads_2[0])|set(tads_2[1])
        intersectsize = len(set(rep1) & set(rep2))
        unionsize = len(set(rep1) | set(rep2))
        return intersectsize / unionsize
    else:
        return 0
'''
for caller in callers:
    print(caller)
    comb = itertools.combinations(li, 2)
    rep_cnt=0
    for i in comb:
        group1=bio1[i[0]]
        group2 = bio1[i[1]]
        for grp1 in group1:
            for grp2 in group2:
                # print(grp1+'--'+grp2)
                ji_tmp=[]
                for chrom in chrs2:
                    # print(caller+'-'+grp1+'-'+grp2+'-'+chrom)

                    current_ji=calculate_JI('../all_TADs/bin/%s/%s_%s.chr%s'%(caller,grp1,caller,chrom),'../all_TADs/bin/%s/%s_%s.chr%s'%(caller,grp2,caller,chrom))
                    # ji_tmp.append(current_ji)
                    JI_c1_bio[rep_cnt][callers.index(caller)]=current_ji
                    rep_cnt+=1
                # JI_c1_bio[rep_cnt][callers.index(caller)]=np.mean(ji_tmp)
                # rep_cnt+=1

    c1_tec_row=c1_tec_rep.shape[0]
    rep_cnt=0
    for c_row in range(0,c1_tec_row):
        ji_tmp = []
        for chrom in chrs2:
            # print(caller+'-'+c1_tec_rep[c_row][0]+'-'+c1_tec_rep[c_row][1]+'-'+chrom)
            current_ji =calculate_JI('../all_TADs/bin/%s/%s_%s.chr%s'%(caller,c1_tec_rep[c_row][0],caller,chrom),'../all_TADs/bin/%s/%s_%s.chr%s'%(caller,c1_tec_rep[c_row][1],caller,chrom) )

            # ji_tmp.append(current_ji)
            JI_c1_tec[rep_cnt][callers.index(caller)] = current_ji
            rep_cnt+=1

        # JI_c1_tec[c_row][callers.index(caller)] = np.mean(ji_tmp)

'''
import seaborn as sn
sn.set(style="ticks", palette="muted", color_codes=True)    

'''
JI_c1_bio=pd.DataFrame(JI_c1_bio,columns=Case_sensitive_callers)
JI_c1_bio['group']='bio'
JI_c1_tec=pd.DataFrame(JI_c1_tec,columns=Case_sensitive_callers)
JI_c1_tec['group']='tec'
JI_c1=pd.concat([JI_c1_bio,JI_c1_tec],ignore_index=True)
JI_c1=pd.melt(JI_c1,id_vars=['group'],value_vars=Case_sensitive_callers,var_name='callers')
sn.boxplot(x='callers',y='value',data=JI_c1,hue='group',ax=ax1,showfliers=False,palette=sn.color_palette("Blues",2))

ax1.set_ylabel('Jaccard Index ',size=13)
ax1.set_xlabel('')
ax1.legend(loc='right', bbox_to_anchor=(1.1, 0.80))
# ax1.tick_params(labelsize=16)
# ax1.legend(prop={'size': 13})
for tick in ax1.get_xticklabels():
    tick.set_rotation(60)

ax1.text(-0.1,1.05, 'A', transform=ax1.transAxes,fontdict = {'size': 22, 'color': 'black'})
'''

largest_span = 0
res_acc_ratio = pd.DataFrame(columns=('callers', 'rep', 'value'))
for caller in callers:
    print(caller)
    repro = pd.read_csv('%s_%s.txt'%(tab,caller), sep='\t', header=None, usecols=[3])
    row_repro = repro.shape[0]
    # print(list(repro[3]).count(1))
    unique_count = repro.nunique()[3]
    largest_span = unique_count if largest_span < unique_count else largest_span
    cnt_list = []
    cnt_sum = []
    for i in range(1, unique_count + 1):
        #current_cnt = repro[3].apply(lambda x: 1 if x == i else 0).sum()
        current_cnt=repro[repro[3]==i].shape[0]
        cnt_list.append((current_cnt / row_repro)*100)
    for i in range(1, unique_count + 1):
        #cnt_sum.append(np.sum(cnt_list[i - 1:unique_count]))
        res_acc_ratio = res_acc_ratio.append([{'callers':caller,'rep':i,'value':np.sum(cnt_list[i - 1:unique_count])}], ignore_index=True)

    #ax3.plot(list(range(1, unique_count + 1)), cnt_sum, label=Case_sensitive_callers[callers.index(caller)],
    #         linewidth=1, color=colors[callers.index(caller)],marker='o', markersize=2)

sn.barplot(x='callers', y='value',hue='rep',data=res_acc_ratio,ax=ax3,palette=sn.color_palette("Blues",4) )

# ax3.set_xticks(ticks)
# ticks = ['≥'+str(x) for x in ticks]
# ax3.set_xticklabels(ticks)
# ax1.legend()
ax3.set_ylabel('''Accumulated
ratio(%)''', size=13)
ax3.set_xlabel('')
ax3.legend(loc='right', bbox_to_anchor=(1.08, 0.5))

# ax1.legend(prop={'size': 11})
# for tick in ax3.get_xticklabels():
#     tick.set_rotation(90)
ax3.text(-0.1,1.05, 'B', transform=ax3.transAxes,fontdict = {'size': 22, 'color': 'black'})
# ax3.legend(loc='right', bbox_to_anchor=(2, 0.45), ncol=2,prop={'size': 12})


# reads_rep1=pd.read_csv('../3_5_data_statistics/res_reads_replicate1.txt', sep='\t')

reads_rep1=None
for caller in callers:
    cur_contact=pd.read_csv('rep_size_contacts/contact_'+caller+'.txt', sep='\t')
    if reads_rep1 is None:
        reads_rep1=cur_contact
    else:
        reads_rep1 = reads_rep1.append(cur_contact, ignore_index=True)


reads_rep1['caller'] = reads_rep1['caller'].map({'Armatus': 'armatus', 'Arrowhead': 'arrowhead',
'CaTCH':'catch','CHAC':'chac', 'CHDF': 'chdf', 'ClusterTAD': 'clustertad', 'deDoc': 'dedoc', 'DI': 'di',
'EAST':'east','GMAP':'gmap','HiCExplorer':'hicexplorer','HiCseg': 'hicseg', 'IC-Finder': 'ic-finder', 
'InsulationScore': 'insulationscore', 'Matryoshka':'matryoshka','MrTADFinder':'mrtadfinder','MSTD': 'mstd',
'OnTAD': 'ontad', 'Spectral':'spectral','SpectralTAD': 'spectraltad', 'TADbit':'tadbit', 'TADBD':'tadbd', 
'TADtree': 'tadtree', 'TopDom': 'topdom'})

reads_rep1 = reads_rep1.sort_values(by=['caller', 'freq'], ascending=[True, True])

reads_rep1['caller']=reads_rep1['caller'].map({'armatus': 'Armatus', 'arrowhead': 'Arrowhead',
'catch':'CaTCH','chac':'Constrained HAC', 'chdf': 'CHDF', 'clustertad': 'ClusterTAD', 'dedoc': 'deDoc', 'di': 'DI',
'east':'EAST','gmap':'GMAP','hicexplorer':'HiCExplorer','hicseg': 'HiCseg', 'ic-finder': 'IC-Finder', 
'insulationscore': 'InsulationScore', 'matryoshka':'Matryoshka','mrtadfinder':'MrTADFinder','mstd': 'MSTD',
'ontad': 'OnTAD', 'spectral':'Spectral','spectraltad': 'SpectralTAD', 'tadbd':'TADBD', 'tadbit':'TADbit',
'tadtree': 'TADtree', 'topdom': 'TopDom'})
# reads_rep1 = reads_rep1[reads_rep1.value < 2000]

#reads_rep1['freq']=reads_rep1['freq'].map({1:'1-2',2:'1-2',3:'3-4',4:'3-4',5:'5-6',6:'5-6',7:'7-8',8:'7-8'})

# reads_rep1['freq']=reads_rep1['freq'].map({1:'1-7',2:'1-7',3:'1-7',4:'1-7',5:'1-7',6:'1-7',7:'1-7',8:'8-14',
#                                            9:'8-14',10:'8-14',11:'8-14',12:'8-14',13:'8-14',14:'8-14',15:'15-21',
#                                            16:'15-21',17:'15-21',18:'15-21',19:'15-21',20:'15-21',21:'15-21',22:'22-29',
#                                            23:'22-29',24:'22-29',25:'22-29',26:'22-29',27:'22-29',28:'22-29',29:'22-29'})

# reads_rep1['freq']=reads_rep1['freq'].map({1:'1-10',2:'1-10',3:'1-10',4:'1-10',5:'1-10',6:'1-10',7:'1-10',8:'1-10',
#                                            9:'1-10',10:'1-10',11:'11-20',12:'11-20',13:'11-20',14:'11-20',15:'11-20',
#                                            16:'11-20',17:'11-20',18:'11-20',19:'11-20',20:'11-20',21:'21-29',22:'21-29',
#                                            23:'21-29',24:'21-29',25:'21-29',26:'21-29',27:'21-29',28:'21-29',29:'21-29'})


# sn.boxplot(x='freq', y='value', data=reads_rep1, hue='caller', showfliers=False,ax=ax4, palette=colors)
sn.boxplot(x='caller', y='value', data=reads_rep1, hue='freq', showfliers=False,ax=ax4,palette=sn.color_palette("Blues",4),linewidth=1)
for tick in ax4.get_xticklabels():
    tick.set_rotation(60)
# ax4.legend_.remove()
ax4.set_ylabel('''TAD average 
contact frequencies''',size=13)
ax4.set_xlabel('')
ax4.legend(loc='right', bbox_to_anchor=(1.08, 0.65))
# ax4.tick_params(labelsize=16)
ax4.text(-0.1,1.05, 'C', transform=ax4.transAxes,fontdict = {'size': 22, 'color': 'black'})

# size_rep1=pd.read_csv('../3_5_data_statistics/res_size_replicate1.txt', sep='\t')

size_rep1=None
for caller in callers:
    cur_size=pd.read_csv('rep_size_contacts/size_'+caller+'.txt', sep='\t')
    if size_rep1 is None:
        size_rep1=cur_size
    else:
        size_rep1 = size_rep1.append(cur_size, ignore_index=True)

size_rep1['caller'] = size_rep1['caller'].map({'Armatus': 'armatus', 'Arrowhead': 'arrowhead',
'CaTCH':'catch','CHAC':'chac', 'CHDF': 'chdf', 'ClusterTAD': 'clustertad', 'deDoc': 'dedoc', 'DI': 'di',
'EAST':'east','GMAP':'gmap','HiCExplorer':'hicexplorer','HiCseg': 'hicseg', 'IC-Finder': 'ic-finder', 
'InsulationScore': 'insulationscore', 'Matryoshka':'matryoshka','MrTADFinder':'mrtadfinder','MSTD': 'mstd',
'OnTAD': 'ontad', 'Spectral':'spectral','SpectralTAD': 'spectraltad', 'TADbit':'tadbit', 'TADBD':'tadbd', 
'TADtree': 'tadtree', 'TopDom': 'topdom'})

size_rep1=size_rep1.sort_values(by=['caller','freq'],ascending=[True,True])

size_rep1['caller'] = size_rep1['caller'].map({'armatus': 'Armatus', 'arrowhead': 'Arrowhead',
'catch':'CaTCH','chac':'Constrained HAC', 'chdf': 'CHDF', 'clustertad': 'ClusterTAD', 'dedoc': 'deDoc', 'di': 'DI',
'east':'EAST','gmap':'GMAP','hicexplorer':'HiCExplorer','hicseg': 'HiCseg', 'ic-finder': 'IC-Finder', 
'insulationscore': 'InsulationScore', 'matryoshka':'Matryoshka','mrtadfinder':'MrTADFinder','mstd': 'MSTD',
'ontad': 'OnTAD', 'spectral':'Spectral','spectraltad': 'SpectralTAD', 'tadbd':'TADBD', 'tadbit':'TADbit',
'tadtree': 'TADtree', 'topdom': 'TopDom'})


#size_rep1['freq']=size_rep1['freq'].map({1:'1-2',2:'1-2',3:'3-4',4:'3-4',5:'5-6',6:'5-6',7:'7-8',8:'7-8'})

# size_rep1['freq']=size_rep1['freq'].map({1:'1-7',2:'1-7',3:'1-7',4:'1-7',5:'1-7',6:'1-7',7:'1-7',8:'8-14',
#                                            9:'8-14',10:'8-14',11:'8-14',12:'8-14',13:'8-14',14:'8-14',15:'15-21',
#                                            16:'15-21',17:'15-21',18:'15-21',19:'15-21',20:'15-21',21:'15-21',22:'22-29',
#                                            23:'22-29',24:'22-29',25:'22-29',26:'22-29',27:'22-29',28:'22-29',29:'22-29'})

# size_rep1['freq']=size_rep1['freq'].map({1:'1-10',2:'1-10',3:'1-10',4:'1-10',5:'1-10',6:'1-10',7:'1-10',8:'1-10',
#                                            9:'1-10',10:'1-10',11:'11-20',12:'11-20',13:'11-20',14:'11-20',15:'11-20',
#                                            16:'11-20',17:'11-20',18:'11-20',19:'11-20',20:'11-20',21:'21-29',22:'21-29',
#                                            23:'21-29',24:'21-29',25:'21-29',26:'21-29',27:'21-29',28:'21-29',29:'21-29'})


# sn.boxplot(x='freq', y='value', data=size_rep1, hue='caller', showfliers=False,ax=ax5, palette=colors)
sn.boxplot(x='caller', y='value', data=size_rep1, hue='freq', showfliers=False,ax=ax5,palette=sn.color_palette("Blues",4),linewidth=1)
for tick in ax5.get_xticklabels():
    tick.set_rotation(60)
# box = ax5.get_position()
# ax5.set_position([box.x0, box.y0, box.width*0.80, box.height])
# ax5.legend(prop={'size': 16})
# ax5.legend_.remove()
ax5.set_ylabel('TAD size',size=13)
ax5.set_xlabel('')
ax5.legend(loc='right', bbox_to_anchor=(1.08, 0.65))
ax5.text(-0.1,1.05, 'D', transform=ax5.transAxes,fontdict = {'size': 22, 'color': 'black'})

# ax5.tick_params(labelsize=16)

# patches = [ mpatches.Patch(color=colors[i], label="{:s}".format(callers[i]) ) for i in range(len(colors)) ]
# fig.legend(loc='right',prop={'size': 14},bbox_to_anchor=(1.12,0.6),ncol=1,frameon=False,shadow=False,handles=patches)


ax2.set_xticks([])
ax3.set_xticks([])
ax4.set_xticks([])


plt.tight_layout()
plt.savefig('Supp_Figure7.jpg',dpi=300,bbox_inches = 'tight')

for caller in callers:
    if(os.path.exists(tab + '_' + caller + '.txt')):
        os.remove(tab + '_' + caller + '.txt')
