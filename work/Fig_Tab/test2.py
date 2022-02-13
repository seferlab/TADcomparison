import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


import seaborn as sn
sn.set(style="ticks", palette="muted", color_codes=True)   

plt.close('all')
fig = plt.subplots()
ax3 = plt.subplot(1,1,1)


callers=['Armatus','Arrowhead','CaTCH','CHAC','CHDF','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','InsulationScore','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADbit','TADtree','TopDom']
Case_sensitive_callers=['Armatus','Arrowhead','CaTCH','Constrained HAC','CHDF','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','InsulationScore','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADbit','TADtree','TopDom']

# colors=['darkslategray','darkcyan','navy','c','cyan','lightgreen','lime','green','black','lightcoral','orange','darkorchid','tomato','chocolate','dodgerblue','darkred']
rainbowcolor=iter(plt.cm.rainbow(np.linspace(0,1,24)))
colors=[]
for i in range(24):
    c=next(rainbowcolor)
    colors.append(c)

import itertools

# tmp1=list(range(3,19))
# tmp1=['HIC%03d'%(x) for x in tmp1]
# tmp2=itertools.combinations(tmp1, 2)  
# tmp3=np.zeros(shape=(0,2))
# for i in tmp2:
#     tmp3=np.row_stack((tmp3,[i[0],i[1]]))

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



chrs=list(range(1,23))
chrs = [str(x) for x in chrs]
chrs.append('X')


dataset=list(range(1,30))
dataset=['HIC%03d'%(x) for x in dataset]

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

'''
#largest_span = 0
res_ratio = pd.DataFrame(columns=('callers', 'rep', 'value'))
for caller in callers:
    repro = pd.read_csv( '%s_%s.txt'%(tab,caller), sep='\t', header=None, usecols=[3])
    row_repro = repro.shape[0]
    # print(list(repro[3]).count(1))
    unique_count = repro.nunique()[3]
    #largest_span = unique_count if largest_span < unique_count else largest_span
    #cnt_list = []
    for i in range(1, unique_count + 1):
        #current_cnt = repro[3].apply(lambda x: 1 if x == i else 0).sum()
        current_cnt=repro[repro[3]==i].shape[0]
        res_ratio = res_ratio.append([{'callers':caller,'rep':i,'value':(current_cnt / row_repro)*100}], ignore_index=True)
        #cnt_list.append((current_cnt / row_repro)*100)
    #ax1.plot(list(range(1, unique_count + 1)), cnt_list, label=Case_sensitive_callers[callers.index(caller)],
    #         linewidth=1, color=colors[callers.index(caller)],marker='o', markersize=2)


print(res_ratio[res_ratio['callers']=='deDoc'])

res_ratio=res_ratio.set_index(['callers','rep'])['value'].unstack()
res_ratio=res_ratio.reindex(index=callers)

res_ratio.plot(kind='bar',stacked = True,ax=ax1,alpha =0.6)
ax1.set_ylim([20,100])
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
    for i in range(1, unique_count + 1, 6):
        #cnt_sum.append(np.sum(cnt_list[i - 1:unique_count]))
        res_acc_ratio = res_acc_ratio.append([{'callers':caller,'rep':i,'value':np.sum(cnt_list[i - 1:unique_count])}], ignore_index=True)

    #ax3.plot(list(range(1, unique_count + 1)), cnt_sum, label=Case_sensitive_callers[callers.index(caller)],
    #         linewidth=1, color=colors[callers.index(caller)],marker='o', markersize=2)

sn.barplot(x='callers', y='value',hue='rep',data=res_acc_ratio,ax=ax3,palette=sn.color_palette("Blues",len(list(range(1, unique_count + 1, 6)))) )

#ticks = list(range(1, largest_span + 1, 1))
#ax3.set_xticks(ticks)
# ticks = ['≥'+str(x) for x in ticks]
# ax3.set_xticklabels(ticks)
# ax1.legend()
ax3.set_ylabel('Accumulated ratio(%)', size=16)
ax3.tick_params(labelsize=16)
#ax3.xaxis.set_major_locator(x_major_locator)
#ax3.set_xlim(0.5,29.5)

# ax1.legend(prop={'size': 11})
for tick in ax3.get_xticklabels():
    tick.set_rotation(60)
ax3.text(-0.1,1.05, 'c', transform=ax3.transAxes,fontdict = {'size': 22, 'color': 'black'})
# ax3.legend(loc='right', bbox_to_anchor=(2, 0.45), ncol=2,prop={'size': 12})








#res_ratio.set_index(['callers','rep'])['value'].unstack().plot.barh(ax = ax1,grid = True,stacked=True,colormap = 'Blues_r')
'''
#ax1.set_xticks(np.arange(1, largest_span + 1, 1))
ax1.set_ylabel('Ratio(%)', size=16)
ax1.tick_params(labelsize=16)

# ax3.legend(loc="upper right",prop={'size': 13})
for tick in ax1.get_xticklabels():
    tick.set_rotation(90)

#ax1.text(-0.1,1.05, 'b', transform=ax1.transAxes,fontdict = {'size': 22, 'color': 'black'})
'''





'''
colors=sn.color_palette("Blues",33)
colors=colors[-29:]
sn.set_palette(colors)
'''


#JI_c1_bio=pd.DataFrame(JI_c1_bio,columns=Case_sensitive_callers)
#JI_c1_bio['group']='bio'
#JI_c1_tec=pd.DataFrame(JI_c1_tec,columns=Case_sensitive_callers)
#JI_c1_tec['group']='tec'
#JI_c1=pd.concat([JI_c1_bio,JI_c1_tec],ignore_index=True)
#JI_c1=pd.melt(JI_c1,id_vars=['group'],value_vars=Case_sensitive_callers,var_name='callers')
#sn.boxplot(x='callers',y='value',data=JI_c1,hue='group',ax=ax1,showfliers=False)

#ax1.set_ylabel('Jaccard Index ',size=16)
#ax1.set_xlabel('')
#ax1.tick_params(labelsize=16)
#ax1.legend(prop={'size': 13})
#for tick in ax1.get_xticklabels():
#    tick.set_rotation(90)

#ax1.text(-0.1,1.05, 'a', transform=ax1.transAxes,fontdict = {'size': 22, 'color': 'black'})


plt.tight_layout()
plt.show()


