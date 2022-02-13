import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
import sys

# callers=['armatus','arrowhead','CHDF','ClusterTAD','deDoc','di','hicseg','IC-Finder','insulation','MSTD','OnTAD','SpectralTAD','tadtree','topdom']

chrs=list(range(1,23))
chrs = [str(x) for x in chrs]
chrs.append('X')

def linear_normalization(data):
    data_shape = data.shape
    data_rows = data_shape[0]
    data_cols = data_shape[1]
    max_x=np.max(data)
    min_x=np.min(data)
    scale_x=max_x-min_x

    for i in range(0, data_rows):
        for j in range(i, data_cols):
            if i==j:
                data[i][j]=(data[i][j]-min_x)/scale_x
            elif i!=j:
                tmp=(data[i][j]-min_x)/scale_x
                data[i][j] = tmp
                data[j][i] = tmp
            

    return data

def reproducibility_frequency(tab,dataset,caller):
    res_size = pd.DataFrame(columns=('caller', 'freq', 'value'))
    res_reads = pd.DataFrame(columns=('caller', 'freq', 'value'))
    TAD_dict = {}
    TAD_size = {}
    TAD_reads = {}
    print(caller)
    for data in dataset:
        print(data)
        for chr in chrs:
            print(chr)
            current_TAD = pd.read_csv('../all_TADs/bin/%s/%s_%s.chr%s'%(caller,data,caller,chr),sep='\t', header=None)
            current_TAD=current_TAD.drop_duplicates()
            current_mat = np.loadtxt('../Rao/%s/%s_50k_KR.chr%s'%(data,data,chr),delimiter='\t')
            current_mat=linear_normalization(current_mat)	
            mat_rows = current_mat.shape[0]
            for index in range(0, current_TAD.shape[0]):
                #start_d = int(current_TAD.iat[index, 1]/40000)
                #right_bound = int(current_TAD.iat[index, 2]/40000)
                start_d = current_TAD.iat[index, 0]
                right_bound = current_TAD.iat[index, 1]
                end_d = right_bound if right_bound < mat_rows else mat_rows
                select_tad = current_mat[start_d:end_d, start_d:end_d]

                if ('chr%s'%(chr),current_TAD.iat[index, 0], current_TAD.iat[index, 1]) in TAD_dict.keys():
                    TAD_dict[('chr%s'%(chr),current_TAD.iat[index, 0], current_TAD.iat[index, 1])] += 1
                    TAD_size[('chr%s'%(chr),current_TAD.iat[index, 0], current_TAD.iat[index, 1])].append((current_TAD.iat[index, 1]-current_TAD.iat[index, 0]+1))
                    TAD_reads[('chr%s'%(chr),current_TAD.iat[index, 0], current_TAD.iat[index, 1])].append(np.mean(np.triu(select_tad)))

                else:
                    TAD_dict[('chr%s'%(chr),current_TAD.iat[index, 0], current_TAD.iat[index, 1])] = 1
                    TAD_size[('chr%s'%(chr),current_TAD.iat[index, 0], current_TAD.iat[index, 1])] = [(current_TAD.iat[index, 1]-current_TAD.iat[index, 0]+1)]
                    TAD_reads[('chr%s'%(chr),current_TAD.iat[index, 0], current_TAD.iat[index, 1])]=[np.mean(np.triu(select_tad))]
    TAD_dict_res = [[x[0], x[1], x[2], TAD_dict[x]] for x in TAD_dict.keys()]
    TAD_dict_res = pd.DataFrame(TAD_dict_res)
    TAD_dict_res = TAD_dict_res.sort_values(by=[0, 3], ascending=True)
    TAD_dict_res.to_csv('rep_size_contacts/%s_%s.txt'%(tab,caller), sep='\t', header=None, index=None)
    tmp_cnt=0
    for key in TAD_dict.keys():
        tmp_cnt+=1
        # sys.stdout.write("\rHandling %d" % (tmp_cnt))
        # sys.stdout.flush()
        for val in TAD_size[key]:
            res_size = res_size.append(pd.DataFrame({'caller': [caller], 'freq': [TAD_dict[key]], 'value': [val]}),ignore_index=True)

    tmp_cnt = 0
    for key in TAD_dict.keys():
        tmp_cnt += 1
        # sys.stdout.write("\nHandling %d" % (tmp_cnt))
        # sys.stdout.flush()
        for val in TAD_reads[key]:
            res_reads = res_reads.append(pd.DataFrame({'caller': [caller], 'freq': [TAD_dict[key]], 'value': [val]}), ignore_index=True)


    res_size.to_csv('rep_size_contacts/size_%s.txt'%(caller), sep='\t', header=True, index=None)
    res_reads.to_csv('rep_size_contacts/contact_%s.txt'%(caller), sep='\t', header=True, index=None)


#reps=list(range(1,30))
#reps=['HIC%03d'%(x) for x in reps]

reps=['HIC003','HIC005','HIC014','HIC020','HIC022','HIC023','HIC025','HIC026']

reproducibility_frequency('replicate1',reps,sys.argv[1])



