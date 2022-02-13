import numpy as np
import pandas as pd
import scipy.stats as st
import math
#import matplotlib.pyplot as plt
import sys


def enrichment_analysis(chrname,chrsize,resolution_kb,TADfile,peak_file):
    res_step = 10
    base_interval=int(1000/res_step) #seperate 1000kb into 10kb intervals
    base_unit=int(1000*res_step) #10kb intervals

    Tads_chr=pd.read_csv(TADfile, sep='\t', header=None)
    # print(Tads_chr.iloc[0])
    unique_boundary=set(Tads_chr[1])|set(Tads_chr[2])
    Peaks=pd.read_csv(peak_file, sep='\t', header=None)
    Peaks_chr= Peaks.iloc[:,1:3].loc[Peaks[0]==chrname]
    Enrich_all = None

    extended_boundary = pd.DataFrame(columns=['left', 'right'])
    # print(len(unique_boundary))
    for boundary in unique_boundary:
        left_bound=boundary-500*1000
        right_bound=boundary+500*1000
        extended_boundary = extended_boundary.append({'left': left_bound, 'right': right_bound}, ignore_index=True)
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
                # print(Enrich_tmp)
                # print(Enrich_tmp.shape)
            if Enrich_all is None:
                Enrich_all=np.array(Vec_tad).reshape(base_interval,1)
            else:
                Enrich_all=np.hstack((Enrich_all,np.array(Vec_tad).reshape(base_interval,1)))
    covered_peak = 0
    for pi in range(0, Peaks_chr.shape[0]):
        left_peak = Peaks_chr.iat[pi, 0]
        right_peak = Peaks_chr.iat[pi, 1]
        boundary_peak = extended_boundary.loc[(extended_boundary['right'] >= left_peak) & (extended_boundary['left'] <= right_peak)]
        if boundary_peak.shape[0] > 0:
            covered_peak += 1
    

    if Enrich_all is None:
        Enrich_all=np.zeros((base_interval,1))
        enrich_line = Enrich_all.mean(axis=1)
        return 0,0,0,0,enrich_line
    else:
        # print(Enrich_all.shape)
        ratio_tmp=Enrich_all[int(base_interval/2-resolution_kb/res_step):int(base_interval/2+resolution_kb/res_step)].sum(axis=0)
        ratio_tmp2 = [1 for i in ratio_tmp.tolist() if i>0 ]
        #boundary_tagged_ratio=sum(ratio_tmp2)/Enrich_all.shape[1]
        
        #ji_element = covered_peak / (Peaks_chr.shape[0] + Enrich_all.shape[1])
        #ji_border = sum(ratio_tmp2) / (Peaks_chr.shape[0] + Enrich_all.shape[1])
        # print(boundary_tagged_ratio)

        enrich_line=Enrich_all.mean(axis=1)
        return covered_peak, sum(ratio_tmp2), (Peaks_chr.shape[0] + Enrich_all.shape[1]), Enrich_all.shape[1], enrich_line
        #return ji_element,ji_border,boundary_tagged_ratio,enrich_line

def extract(element,peaks_bedfile):
    res_step = 10
    base_interval = int(1000 / res_step)  # seperate 1000kb into 10kb intervals
    callers=['Armatus','Arrowhead','CaTCH','CHAC','CHDF','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','InsulationScore','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADbit','TADtree','TopDom']

    dataset=list(range(1,30))
    dataset=['HIC%03d'%(x) for x in dataset]
    chrs=list(range(1,23))
    chrs = [str(x) for x in chrs]
    chrs.append('X')
    chrs_length=[249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,
                 141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,
                 81195210,78077248,59128983,63025520,48129895,51304566,155270560]

    overall_res=np.zeros((len(callers),6,len(dataset)))
    comparing_items=['Peak maximum','Boundary tagged ratio','Fold change','JI_element','JI_border','P-value ratio']

    for caller in callers:
        print(caller)
        res_caller=np.zeros((6,len(dataset)))
        for data in dataset:
            plot_lines=np.zeros((100,23))
            covered_element = []
            covered_boundary = []
            boundary_element_assemble = []
            boundary_assemble = []
            for chr in chrs:
                curr_covered_element,curr_covered_boundary,curr_BEA,curr_boundary,plot_lines[:,chrs.index(chr)]=enrichment_analysis('chr'+chr, chrs_length[chrs.index(chr)], 50,
                                 '../all_TADs/loci/%s/%s_%s.chr%s'%(caller,data,caller,chr),peaks_bedfile)
                covered_element.append(curr_covered_element)
                covered_boundary.append(curr_covered_boundary)
                boundary_element_assemble.append(curr_BEA)
                boundary_assemble.append(curr_boundary)
                
            mean_peak = plot_lines.mean(axis=1)
            peak_maximum = np.max(mean_peak[int(base_interval / 2 - base_interval * 0.01):int(base_interval / 2 + base_interval * 0.01)])
            # print(peak_maximum)
            peak_val = np.mean(mean_peak[int(base_interval / 2 - base_interval * 0.01):int(base_interval / 2 + base_interval * 0.01)])
            bg = mean_peak[0:int(base_interval / 2 - base_interval * 0.2)].tolist() + mean_peak[int(base_interval / 2 + base_interval * 0.2 - 1):base_interval].tolist()
            bg_val = np.mean(bg)
            fc = math.log2(peak_val / bg_val)
            # print(fc)
            pval_peak = 1-st.norm.cdf(peak_val, loc=bg_val, scale=np.std(bg))
            # print(pval_peak)
            boundary_tagged_ratio=np.sum(covered_boundary)/np.sum(boundary_assemble)
            ji_element=np.sum(covered_element)/np.sum(boundary_element_assemble)
            ji_border=np.sum(covered_boundary)/np.sum(boundary_element_assemble)
            res_caller[:,dataset.index(data)]=[peak_maximum,boundary_tagged_ratio,fc,ji_element,ji_border,pval_peak]


        overall_res[callers.index(caller),:,:]=res_caller

    for i in range(0,5):
        data_to_plot=overall_res[:, i, :].T
        # plt.close('all')
        # fig, ax1 = plt.subplots()
        # ax1.boxplot(data_to_plot)
        # ax1.set_ylabel(comparing_items[i])
        # ax1.set_xticklabels(Case_sensitive_callers)
        # for tick in ax1.get_xticklabels():
        #     tick.set_rotation(90)
        # plt.tight_layout()
        # plt.savefig('./out/'+element+' '+comparing_items[i]+'.pdf')
        np.savetxt('output/%s_%s.txt'%(element,comparing_items[i]), data_to_plot) 

    data_to_plot = overall_res[:, 5, :]
    data_to_plot[data_to_plot>0.05]=1
    data_to_plot[data_to_plot <= 0.05] = 0
    data_to_plot = data_to_plot.mean(axis=1)
    data_to_plot=1-data_to_plot
    # plt.close('all')
    # fig, ax1 = plt.subplots()
    # ax1.plot(Case_sensitive_callers, data_to_plot, linewidth=2, color='b', marker='o', markersize=4)
    # ax1.set_ylabel(comparing_items[3])
    # for tick in ax1.get_xticklabels():
    #     tick.set_rotation(90)
    # plt.tight_layout()
    # plt.savefig('./out/' + element + ' ' + comparing_items[4] + '.pdf')
    np.savetxt('output/%s_%s.txt'%(element,comparing_items[5]), data_to_plot) 


#extract('CTCF','marks/wgEncodeAwgTfbsSydhGm12878Ctcfsc15914c20UniPk.narrowPeak')
#extract('RAD21','marks/wgEncodeAwgTfbsSydhGm12878Rad21IggrabUniPk.narrowPeak')
#extract('SMC3','marks/wgEncodeAwgTfbsSydhGm12878Smc3ab9263IggmusUniPk.narrowPeak')
#extract('H3K4me3','marks/wgEncodeUwHistoneGm12878H3k4me3StdPkRep1.narrowPeak')
#extract('H3K36me3','marks/wgEncodeUwHistoneGm12878H3k36me3StdPkRep1.narrowPeak')
#extract('PolII','marks/wgEncodeOpenChromChipGm12878Pol2Pk.narrowPeak')
#extract('Housekeeping genes','marks/HK.txt')
#extract('SINE','marks/SINE.bed')

extract(sys.argv[1],sys.argv[2])

