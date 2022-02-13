import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sn

sn.set(style="ticks", palette="muted", color_codes=True) 


def linear_normalization(data):
    data_shape = data.shape
    data_rows = data_shape[0]
    data_cols = data_shape[1]
    max_x = np.max(data)
    min_x = np.min(data)
    scale_x = max_x - min_x

    for i in range(0, data_rows):
        for j in range(i, data_cols):
            if i == j:
                data[i][j] = (data[i][j] - min_x) / scale_x
            elif i != j:
                tmp = (data[i][j] - min_x) / scale_x
                data[i][j] = tmp
                data[j][i] = tmp

    return data


def nonlinearity_normalization(data, max_value):
    data_shape = data.shape
    data_rows = data_shape[0]
    data_cols = data_shape[1]
    log_max = math.log(max_value + 1, 10)

    for i in range(0, data_rows):
        for j in range(i, data_cols):
            if i == j:
                data[i][j] = math.log((data[i][j] + 1), 10) / log_max
            elif i != j:
                tmp = math.log((data[i][j] + 1), 10) / log_max
                data[i][j] = tmp
                data[j][i] = tmp
    return data


def plot_mat_comparison(mat, low_bound, high_bound, current_ax, TAD_file, tool, ytick_flag):
    tadset = np.loadtxt(TAD_file, usecols=(0,1), dtype=int)
    tadset = tadset[tadset[:, 0].argsort()]
    if low_bound < 0 or high_bound > mat.shape[0] or high_bound <= low_bound:
        print("bad parameter.")
    else:
        current_ax.set_xlim([low_bound, high_bound])
        current_ax.set_ylim([low_bound, high_bound])
        current_ax.invert_yaxis()
        current_ax.xaxis.tick_top()
        current_ax.set_xlabel(tool)
        for tick in current_ax.get_xticklabels():
            tick.set_rotation(90)
        if ytick_flag == False:
            current_ax.set_yticks([])
        # current_ax.spines['top'].set_visible(False)
        # current_ax.spines['right'].set_visible(False)
        # current_ax.spines['bottom'].set_visible(False)
        # current_ax.spines['left'].set_visible(False)
        max_value_mat = np.max(mat[low_bound:high_bound, low_bound:high_bound])
        ax_show = current_ax.imshow(mat, cmap=plt.cm.Blues, vmax=max_value_mat, vmin=0)

        for i in range(0, len(tadset)):
            tad_start = tadset[i][0] 
            tad_end = tadset[i][1]
            rect = patches.Rectangle((tad_start, tad_start), tad_end - tad_start, tad_end - tad_start, linewidth=1,
                                     edgecolor='r', facecolor='none')
            current_ax.add_patch(rect)
        # if show_colorbar==True:
        #     position = fig.add_axes([0.93, 0.34, 0.005, 0.5])
        #     fig.colorbar(ax_show, cax=position, orientation='vertical')  # 方向
    return ax_show


inputfile = '../Rao/HIC001/HIC001_50k_KR.chr1'
inputdata = np.loadtxt(inputfile)
inputdata = nonlinearity_normalization(inputdata, np.max(inputdata))
# inputdata=linear_normalization(inputdata)

# fig = plt.gcf()
# fig.set_size_inches(7.0/3,7.0/3) #dpi = 300, output = 700*700 pixels

fig = plt.figure(figsize=(13, 6.6))
ax = fig.subplots(3,8)

plot_mat_comparison(inputdata, 1800, 1900, ax[0][0], '../all_TADs/bin/%s/HIC001_%s.chr1'%('Armatus','Armatus'), 'Armatus', True)
plot_mat_comparison(inputdata, 1800, 1900, ax[0][1], '../all_TADs/bin/%s/HIC001_%s.chr1'%('Arrowhead','Arrowhead'), 'Arrowhead', False)
plot_mat_comparison(inputdata, 1800, 1900, ax[0][2], '../all_TADs/bin/%s/HIC001_%s.chr1'%('CaTCH','CaTCH'), 'CaTCH', False)
plot_mat_comparison(inputdata, 1800, 1900, ax[0][3], '../all_TADs/bin/%s/HIC001_%s.chr1'%('CHAC','CHAC'), 'Constrained HAC', False)
plot_mat_comparison(inputdata, 1800, 1900, ax[0][4], '../all_TADs/bin/%s/HIC001_%s.chr1'%('CHDF','CHDF'), 'CHDF', False)
plot_mat_comparison(inputdata, 1800, 1900, ax[0][5], '../all_TADs/bin/%s/HIC001_%s.chr1'%('ClusterTAD','ClusterTAD'), 'ClusterTAD', False)
plot_mat_comparison(inputdata, 1800, 1900, ax[0][6], '../all_TADs/bin/%s/HIC001_%s.chr1'%('deDoc','deDoc'), 'deDoc', False)
plot_mat_comparison(inputdata, 1800, 1900, ax[0][7], '../all_TADs/bin/%s/HIC001_%s.chr1'%('DI','DI'), 'DI', False)
plot_mat_comparison(inputdata, 1800, 1900, ax[1][0], '../all_TADs/bin/%s/HIC001_%s.chr1'%('EAST','EAST'), 'EAST', True)
plot_mat_comparison(inputdata, 1800, 1900, ax[1][1], '../all_TADs/bin/%s/HIC001_%s.chr1'%('GMAP','GMAP'), 'GMAP', False)
plot_mat_comparison(inputdata, 1800, 1900, ax[1][2], '../all_TADs/bin/%s/HIC001_%s.chr1'%('HiCExplorer','HiCExplorer'), 'HiCExplorer', False)
plot_mat_comparison(inputdata, 1800, 1900, ax[1][3], '../all_TADs/bin/%s/HIC001_%s.chr1'%('HiCseg','HiCseg'), 'HiCseg', False)
plot_mat_comparison(inputdata,1800,1900,ax[1][4],'../all_TADs/bin/%s/HIC001_%s.chr1'%('IC-Finder','IC-Finder'),'IC-Finder',False)
plot_mat_comparison(inputdata,1800,1900,ax[1][5],'../all_TADs/bin/%s/HIC001_%s.chr1'%('InsulationScore','InsulationScore'),'InsulationScore',False)
plot_mat_comparison(inputdata,1800,1900,ax[1][6],'../all_TADs/bin/%s/HIC001_%s.chr1'%('Matryoshka','Matryoshka'),'Matryoshka',False)
plot_mat_comparison(inputdata,1800,1900,ax[1][7],'../all_TADs/bin/%s/HIC001_%s.chr1'%('MrTADFinder','MrTADFinder'),'MrTADFinder',False)
plot_mat_comparison(inputdata,1800,1900,ax[2][0],'../all_TADs/bin/%s/HIC001_%s.chr1'%('MSTD','MSTD'),'MSTD',True)
plot_mat_comparison(inputdata,1800,1900,ax[2][1],'../all_TADs/bin/%s/HIC001_%s.chr1'%('OnTAD','OnTAD'),'OnTAD',False)
plot_mat_comparison(inputdata,1800,1900,ax[2][2],'../all_TADs/bin/%s/HIC001_%s.chr1'%('Spectral','Spectral'),'Spectral',False)
plot_mat_comparison(inputdata,1800,1900,ax[2][3],'../all_TADs/bin/%s/HIC001_%s.chr1'%('SpectralTAD','SpectralTAD'),'SpectralTAD',False)
plot_mat_comparison(inputdata,1800,1900,ax[2][4],'../all_TADs/bin/%s/HIC001_%s.chr1'%('TADBD','TADBD'),'TADBD',False)
plot_mat_comparison(inputdata,1800,1900,ax[2][5],'../all_TADs/bin/%s/HIC001_%s.chr1'%('TADbit','TADbit'),'TADbit',False)
plot_mat_comparison(inputdata,1800,1900,ax[2][6],'../all_TADs/bin/%s/HIC001_%s.chr1'%('TADtree','TADtree'),'TADtree',False)
im=plot_mat_comparison(inputdata,1800,1900,ax[2][7],'../all_TADs/bin/%s/HIC001_%s.chr1'%('TopDom','TopDom'),'TopDom',False)

# fig.set_tight_layout(True)
fig.subplots_adjust(left=0.08)
# 左下宽高
# ([0.72, 0.24, 0.01, 0.5])
cbar_ax = fig.add_axes([0.92, 0.24, 0.01, 0.5])
fig.colorbar(im, cax=cbar_ax)
# plt.show()

# plt.tight_layout()#调整整体空白
# plt.savefig('out_1800_1900.pdf')
plt.savefig('Supp_Figure3.jpg',dpi=300, bbox_inches='tight')
# 3000,3100
# 1800,1900
