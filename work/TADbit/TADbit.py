"""
24 Oct 2012
"""
from __future__ import print_function

from os                           import path, listdir
from pytadbit.parsers.hic_parser  import read_matrix
from pytadbit.tadbit_py           import _tadbit_wrapper
from math                         import isnan, sqrt
from scipy.sparse.csr             import csr_matrix
from scipy.stats                  import mannwhitneyu
import numpy as np
def tadbit(x, remove=None, n_cpus=1, verbose=True,
           max_tad_size="max", no_heuristic=0, use_topdom=False, topdom_window=5, **kwargs):
    """
    The TADbit algorithm works on raw chromosome interaction count data.
    The normalization is neither necessary nor recommended,
    since the data is assumed to be discrete counts.

    TADbit is a breakpoint detection algorithm that returns the optimal
    segmentation of the chromosome under BIC-penalized likelihood. The
    model assumes that counts have a Poisson distribution and that the
    expected value of the counts decreases like a power-law with the
    linear distance on the chromosome. This expected value of the counts
    at position (i,j) is corrected by the counts at diagonal positions
    (i,i) and (j,j). This normalizes for different restriction enzyme
    site densities and 'mappability' of the reads in case a bin contains
    repeated regions.

    :param x: a square matrix of interaction counts in the HI-C data or a list
       of such matrices for replicated experiments. The counts must be evenly
       sampled and not normalized. x might be either a list of list, a path to
       a file or a file handler
    :argument 'visibility' norm: kind of normalization to use. Choose between
       'visibility' of 'Imakaev'
    :argument None remove: a python list of lists of booleans mapping positively
       columns to remove (if None only columns with a 0 in the diagonal will be
       removed)
    :param 1 n_cpus: The number of CPUs to allocate to TADbit. If
       n_cpus='max' the total number of CPUs will be used
    :param auto max_tad_size: an integer defining maximum size of TAD. Default
       (auto or max) defines it as the number of rows/columns
    :param False no_heuristic: whether to use or not some heuristics
    :param False use_topdom: whether to use TopDom algorithm to find tads or not (http://www.ncbi.nlm.nih.gov/pubmed/26704975, http://zhoulab.usc.edu/TopDom/)
    :param 5 topdom_window: the window size for topdom algorithm
    :param False get_weights: either to return the weights corresponding to the
       Hi-C count (weights are a normalization dependent of the count of each
       columns)

    :returns: the :py:func:`list` of topologically associated domains'
       boundaries, and the corresponding list associated log likelihoods.
       If no weights are given, it may also return calculated weights.
    """
    nums = [hic_data for hic_data in read_matrix(x, one=False)]

    if not use_topdom:
        size = len(nums[0])
        nums = [num.get_as_tuple() for num in nums]
        if not remove:
            # if not given just remove columns with zero in diagonal
            remove = tuple([0 if nums[0][i*size+i] else 1 for i in range(size)])
        n_cpus = n_cpus if n_cpus != 'max' else 0
        max_tad_size = size if max_tad_size in ["max", "auto"] else max_tad_size
        _, nbks, passages, _, _, bkpts = \
           _tadbit_wrapper(nums,             # list of lists of Hi-C data
                           remove,           # list of columns marking filtered
                           size,             # size of one row/column
                           len(nums),        # number of matrices
                           n_cpus,           # number of threads
                           int(verbose),     # verbose 0/1
                           max_tad_size,     # max_tad_size
                           kwargs.get('ntads', -1) + 1,
                           int(no_heuristic),# heuristic 0/1
                           )

        breaks = [i for i in range(size) if bkpts[i + nbks * size] == 1]
        scores = [p for p in passages if p > 0]

        result = {'start': [], 'end'  : [], 'score': []}
        for brk in range(len(breaks)+1):
            result['start'].append((breaks[brk-1] + 1) if brk > 0 else 0)
            result['end'  ].append(breaks[brk] if brk < len(breaks) else size - 1)
            result['score'].append(scores[brk] if brk < len(breaks) else None)
    else:
        result = {'start': [], 'end'  : [], 'score': [], 'tag': []}

        ret = TopDom(nums[0],window_size=topdom_window)


        for key in sorted(ret):
            result['tag'].append(ret[key]['tag'])
            result['start'].append(ret[key]['start'])
            result['end'].append(ret[key]['end'])
            if ret[key]['tag'] == 'domain':
                result['score'].append(ret[key]['score'])
            else:
                result['score'].append(0)

        max_score = max(result['score'])
        for i in range(len(result['score'])):
            result['score'][i] = 1-int((result['score'][i]/max_score)*10)

    return result
import numpy as np
import sys
#txt=np.loadtxt('test_KR.chr21')
#txt[np.isnan(txt)]=0
#np.savetxt('test.test',txt,fmt='%g',delimiter='\t')
#res=tadbit('test.test', remove=None, n_cpus=2)
res=tadbit(sys.argv[1], remove=None, n_cpus=70)
TADs=np.zeros((len(res['start']),2),dtype=np.int)
for i in range(len(res['start'])):
    TADs[i,0]=res['start'][i]+1
    TADs[i,1]=res['end'][i]+1
np.savetxt(sys.argv[2],TADs,delimiter='\t',fmt='%g')
#print(res)

