import numpy as np
import pandas as pd

def extract(cate):
    Case_sensitive_callers=['Armatus','Arrowhead','CaTCH','CHAC','CHDF','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','InsulationScore','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADbit','TADtree','TopDom']
    regulators=['CTCF','Housekeeping_genes','RAD21','TSS']
    BTR=None
    for regulate in regulators:
        if cate=='Peak maximum' and regulate=='H3K9me3':
            tmp_btr = np.loadtxt('./output/' + regulate + ' Valley minimum.txt')
        else:
            tmp_btr=np.loadtxt('./output/'+regulate+'_'+cate+'.txt')
        tmp_btr = pd.DataFrame(tmp_btr, columns=Case_sensitive_callers)
        tmp_btr['group'] = regulate
        if BTR is None:
            BTR=tmp_btr
        else:
            BTR = pd.concat([BTR, tmp_btr], ignore_index=True)

    # print(BTR)
    BTR=BTR.groupby('group').mean()
    # print(BTR.groupby('group').mean())
    BTR =BTR.round(3)
    # BTR.to_csv('output/3_4_'+cate+'.txt',sep='\t')
    BTR.to_csv('output/1_' + cate + '.csv', sep=',')



def extract2(cate):
    Case_sensitive_callers=['Armatus','Arrowhead','CaTCH','CHAC','CHDF','ClusterTAD','deDoc','DI','EAST','GMAP','HiCExplorer','HiCseg','IC-Finder','InsulationScore','Matryoshka','MrTADFinder','MSTD','OnTAD','Spectral','SpectralTAD','TADBD','TADbit','TADtree','TopDom']
    regulators=['CTCF','Housekeeping_genes','RAD21','TSS']
    BTR=None
    for regulate in regulators:
        tmp_btr=np.loadtxt('./output/'+regulate+'_'+cate+'.txt').reshape(1,len(Case_sensitive_callers))
        # print(tmp_btr)
        tmp_btr = pd.DataFrame(tmp_btr, columns=Case_sensitive_callers)
        tmp_btr['group'] = regulate
        if BTR is None:
            BTR=tmp_btr
        else:
            BTR = pd.concat([BTR, tmp_btr], ignore_index=True)

    # print(BTR)
    BTR = BTR.groupby('group').mean()
    # print(BTR.groupby('group').mean())
    BTR = BTR.round(3)
    # BTR.to_csv('output/3_4_' + cate + '.txt', sep='\t')
    BTR.to_csv('output/1_' + cate + '.csv', sep=',')
    



extract('Boundary tagged ratio')
extract('Fold change')
extract('Peak maximum') #Average peak metric
extract('JI_element')
extract('JI_border')
extract2('P-value ratio')
