import numpy as np
import pandas as pd
import scipy.stats as st
import math

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def TSS_enrichment(chrname,chrsize,TADfile,peak_file):
    res_step = 10
    base_interval=int(1000/res_step) #seperate 1000kb into 10kb intervals
    base_unit=int(1000*res_step) #10kb intervals

    Tads_chr=pd.read_csv(TADfile, sep='\t', header=None)
    # print(Tads_chr.iloc[0])
    unique_boundary=set(Tads_chr[1])|set(Tads_chr[2])
    Peaks = pd.read_csv(peak_file, sep='\t', header=None)
    Peaks_chr = Peaks.iloc[:, 1:3].loc[Peaks[1] == chrname]
    Peaks_chr = Peaks_chr[2].tolist()
    Enrich_all = None

    for boundary in unique_boundary:
        left_bound=boundary-500*1000
        right_bound=boundary+500*1000
        if left_bound > 0 and right_bound < chrsize:
            Peaks_tad = [i for i in Peaks_chr if i >= left_bound and i <= right_bound]
            # Peaks_tad = Peaks_chr.iloc[:, :].loc[(Peaks_chr[3] >= left_bound) & (Peaks_chr[2] <= right_bound)]
            Vec_tad= np.zeros((base_interval, 1), dtype=int).flatten().tolist()
            for pi in range(0, len(Peaks_tad)):
                Enrich_current = np.zeros((base_interval, 1), dtype=int).flatten().tolist()
                current_peak = Peaks_tad[pi]
                if current_peak==right_bound:
                    Enrich_current[99] += 1
                elif current_peak > left_bound:
                    Enrich_current[int((current_peak - left_bound) / base_unit)] += 1
                # if current_peak > left_bound:
                #     Enrich_current[int((current_peak - left_bound) / base_unit)] += 1
                Enrich_tmp = [1 if i > 0 else 0 for i in Enrich_current]
                for i in range(0, base_interval):
                    Vec_tad[i] = Vec_tad[i] + Enrich_tmp[i]
                # print(Enrich_tmp)
                # print(Enrich_tmp.shape)
            if Enrich_all is None:
                Enrich_all=np.array(Vec_tad).reshape(base_interval,1)
            else:
                Enrich_all=np.hstack((Enrich_all,np.array(Vec_tad).reshape(base_interval,1)))

    # enrich_line=Enrich_all.mean(axis=1)
    # return enrich_line
    if Enrich_all is None:
        return [0]*100
    else:
        enrich_line=Enrich_all.mean(axis=1)
        return enrich_line

def cohesin_enrichment(chrname,chrsize,TADfile,peak_file):
    res_step = 10
    base_interval=int(1000/res_step) #seperate 1000kb into 10kb intervals
    base_unit=int(1000*res_step) #10kb intervals

    Tads_chr=pd.read_csv(TADfile, sep='\t', header=None)
    # print(Tads_chr.iloc[0])
    unique_boundary=set(Tads_chr[1])|set(Tads_chr[2])
    Peaks=pd.read_csv(peak_file, sep='\t', header=None)
    Peaks_chr= Peaks.iloc[:,1:3].loc[Peaks[0]==chrname]
    Enrich_all = None

    # print(len(unique_boundary))
    for boundary in unique_boundary:
        left_bound=boundary-500*1000
        right_bound=boundary+500*1000
        if left_bound > 0 and right_bound < chrsize:
            Peaks_tad=Peaks_chr.iloc[:,:].loc[(Peaks_chr[2]>=left_bound) & (Peaks_chr[1]<=right_bound)]
            Vec_tad= np.zeros((base_interval, 1), dtype=int).flatten().tolist()
            for pi in range(0,Peaks_tad.shape[0]):
                Enrich_current = np.zeros((base_interval, 1), dtype=int).flatten().tolist()
                current_peak = Peaks_tad.iat[pi, 0]
                start_peak = current_peak if current_peak > left_bound else left_bound
                current_peak = Peaks_tad.iat[pi, 1]
                end_peak = current_peak if current_peak < right_bound else right_bound - 1
                start_peak = int((start_peak - left_bound) / base_unit)
                end_peak = int((end_peak - left_bound) / base_unit)
                for i in range(start_peak, end_peak + 1):
                    Enrich_current[i] += 1

                Enrich_tmp = [1 if i>0 else 0 for i in Enrich_current]
                for i in range(0, base_interval):
                    Vec_tad[i] = Vec_tad[i] + Enrich_tmp[i]

            if Enrich_all is None:
                Enrich_all=np.array(Vec_tad).reshape(base_interval,1)
            else:
                Enrich_all=np.hstack((Enrich_all,np.array(Vec_tad).reshape(base_interval,1)))
    if Enrich_all is None:
        return [0]*100
    else:
        enrich_line=Enrich_all.mean(axis=1)
        return enrich_line
   
chrs=list(range(1,23))
chrs = [str(x) for x in chrs]
chrs.append('X')
chrs_length=[249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,
             141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,
             81195210,78077248,59128983,63025520,48129895,51304566,155270560]

dataset='HIC001'
callers=['Armatus','Arrowhead','CaTCH','CHAC','CHDF','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','InsulationScore','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADbit','TADtree','TopDom']
Case_sensitive_callers=['Armatus','Arrowhead','CaTCH','''Constrained
HAC''','CHDF','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','''Insulation
Score''','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADbit','TADtree','TopDom']

# colors=['darkslategray','darkcyan','navy','c','cyan','lightgreen','lime','green','black','lightcoral','orange','darkorchid','tomato','chocolate','dodgerblue','darkred']
rainbowcolor=iter(plt.cm.rainbow(np.linspace(0,1,24)))
colors=[]
for i in range(24):
    c=next(rainbowcolor)
    colors.append(c)


regulators=['CTCF','H3k9me3']
regu_files=['marks/wgEncodeAwgTfbsSydhGm12878Ctcfsc15914c20UniPk.narrowPeak',
'marks/wgEncodeBroadHistoneGm12878H3k9me3StdPk.broadPeak']


# fig = plt.figure(figsize=(25,10))
# # ax = fig.subplots(2,5, projection='3d')
# axes=[]
# for i in range(1,11):
#     ax1 = fig.add_subplot(2,5,i, projection='3d')
#     axes.append(ax1)
fig = plt.figure(figsize=(13, 11))
axes = fig.subplots(6,8)


sss='abcdefghij'
color=(0.41708573625528644, 0.6806305267204922, 0.8382314494425221)

handle_legend=[]

for re in range(len(regulators)):
    for caller in range(len(callers)):
        print(caller)
        plot_lines=np.zeros((100,23))
        for chr in chrs:
            plot_lines[:,chrs.index(chr)]=cohesin_enrichment('chr'+chr, chrs_length[chrs.index(chr)], '../all_TADs/loci/%s/HIC001_%s.chr%s'%(callers[caller],callers[caller],chr), regu_files[re])

        # print(plot_lines.mean(axis=1).shape)

        l1, = axes[re*3+int(caller/8)][int(caller%8)].plot(range(0, 100), plot_lines.mean(axis=1), linewidth=1, color=color, marker='o',markersize=0)

        # if re==9:
        #     handle_legend.append(l1)

        # ax1.plot(range(0,100), plot_lines.mean(axis=1), label=regulators[re], linewidth=1, color=colors[re], marker='o', markersize=0)

        axes[re*3+int(caller/8)][int(caller%8)].set_xticks(list(range(0,101,50)))
        # axes[re].set_xticklabels(['-0.5Mb', '-0.25Mb', '0', '+0.25Mb', '+0.5Mb'])
        axes[re*3+int(caller/8)][int(caller%8)].set_yticklabels([])
        axes[re*3+int(caller/8)][int(caller%8)].set_xticklabels(['-500kb', '0', '+500kb'])
        # axes[re*3+int(caller/8)][int(caller%8)].tick_params(labelsize=16)
        axes[re*3+int(caller/8)][int(caller%8)].set_title(callers[caller])
        # ax1.set_ylabel('Average number of elements per 10 kb',size=16)
        # axes[re].set_title(sss[re:re+1],fontdict={'color': 'k', 'fontsize': 22},loc='left')
        # axes[re].text(-0.1, 1.05, sss[re:re+1], transform=axes[re].transAxes, fontdict={'size': 22, 'color': 'black'})

axes[0][0].text(-0.4, 1.1, 'A', transform=axes[0][0].transAxes, fontdict={'size': 22, 'color': 'black'})
axes[3][0].text(-0.4, 1.1, 'B', transform=axes[3][0].transAxes, fontdict={'size': 22, 'color': 'black'})

# fig.legend(handles=handle_legend,loc='right', bbox_to_anchor=(1.11, 0.48),fontsize=18)

plt.tight_layout()
plt.savefig('Supp_Figure4.jpg',dpi=300,bbox_inches = 'tight')
# plt.show()

