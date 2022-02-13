import numpy as np
from timeit import time
from numba import guvectorize,float64,int64
import os
from scipy import sparse
from detect_peaks import detect_peaks
import importlib
pnd = importlib.find_loader('pandas')
PANDAS_INSTALLED = pnd is not None
if PANDAS_INSTALLED:
    import pandas as pd
import sys

RESOLUTION = int(sys.argv[3])
if RESOLUTION < 500000:
    W = 20 # insulation window size
else:
    W = 10
maxL = 2*int(np.round(3200000/RESOLUTION)) + 1 # maximum length of TAD allowed
Nfactor = 0.35 # normalization factor: larger values lead to smaller TADs



st = time.time()
if PANDAS_INSTALLED:
    chr1 = pd.read_csv(sys.argv[1],sep='\t',header=None)
    chr1 = chr1.values
else:
    chr1 = np.genfromtxt(sys.argv[1])
print('time to read the input data :',time.time()-st)
N = chr1.shape[0]

# compute the integral image (We ignore the values on the diagonal)
st = time.time()
intgMat = np.zeros([2*maxL,N],dtype=np.float64)
I = np.arange(N-1) # delta = 1                
intgMat[1,I+1] = chr1[I,I+1] + intgMat[0,I+1] + intgMat[0,I] 
for delta in range(2,2*maxL):
    I = np.arange(N-delta)
    intgMat[delta,I+delta] = chr1[I,I+delta] + intgMat[delta-1,I+delta] + intgMat[delta-1,I+delta-1] - intgMat[delta-2,I+delta-1] 
print('time to compute the integral image:',time.time()-st)
#print('********************** TAD DETECTION *****************************')
# TAD Detection
#@guvectorize([(float64[:,:], int64[:], int64[:], float64[:])], '(m,p),(),()->()',target='parallel',nopython=False)
def score(intgMAT,i,l):
    i_indent = maxL
    j_indent = maxL
    wScore = 0
    pixel = i[0]
    if(l[0]<=maxL or l[0]<5):
        if l[0] % 2 == 0:
            
            w = (l[0])/2
            pixel = pixel + w
            w2 = np.math.ceil(w/5)

            A = [int(pixel-w-l[0]+i_indent),int(pixel-w+j_indent)]
            B = [int(pixel-w-l[0]+i_indent),int(pixel+w+j_indent)]
            D = [int(pixel-w+i_indent),int(pixel+w+j_indent)]
            E = [int(pixel-w+i_indent),int(pixel+w+l[0]+j_indent)]
            F = [int(pixel+w+i_indent),int(pixel+w+l[0]+j_indent)]
        else:
            w = (l[0]-1)/2
            pixel = pixel + w
            w2 = np.math.ceil(w/5)

            A = [int(pixel-w-l[0]-1+i_indent),int(pixel-w-1+j_indent)]
            B = [int(pixel-w-l[0]-1+i_indent),int(pixel+w+j_indent)]
            D = [int(pixel-w-1+i_indent),int(pixel+w+j_indent)]
            E = [int(pixel-w-1+i_indent),int(pixel+w+l[0]+j_indent)]
            F = [int(pixel+w+i_indent),int(pixel+w+l[0]+j_indent)]

        wScore = intgMAT[D[1]-D[0]+i_indent, D[1]]
        #res[0] = wScore/np.power(l[0],Nfactor)
        return wScore/np.power(l[0],Nfactor)
    else:
        #res[0]=0
        return 0

@guvectorize([(float64[:,:], int64[:], int64[:], float64[:])], '(m,p),(),(o)->(o)',target='parallel')
def det_score(intgMAT,i,w,res):
    i_indent = maxL
    j_indent = maxL
    
    B = [int(i[0]-w[0]+i_indent),int(i[0]+j_indent)]
    C = [int(i[0]-w[0]+i_indent),int(i[0]+w[0]+j_indent)]
    F = [int(i[0]+i_indent),int(i[0]+w[0]+j_indent)]

    a = intgMAT[B[1]-B[0]+i_indent, B[1]] #left
    f = intgMAT[F[1]-F[0]+i_indent, F[1]] #right
    res[0] = intgMAT[C[1]-C[0] +i_indent, C[1]] - a - f #center
    res[1] = (a - res[0])
    res[2] = (f - res[0])

@guvectorize([(float64[:,:],int64[:],int64[:],float64[:], int64[:], int64[:],int64[:],int64[:],int64[:], float64[:])], '(u,v),(s),(s),(r),(p),(q),(p),(),(n)->(n)',target='parallel')
def parDP(scores,i2I,j2J,T,startPeaks,endPeaks,prev_end,j,k,res):
    for counter in range(k.shape[0]):
        res[counter] = T[prev_end[k[counter]]] + scores[i2I[startPeaks[k[counter]]], j2J[endPeaks[j[0]]]]#scores[startPeaks[k[counter]], endPeaks[j[0]]]
       
def findTADs(intgMAT):
    st = time.time()
    intgMAT2 = np.lib.pad(intgMAT, ((maxL, maxL), (maxL, maxL)), 'constant', constant_values=0)       
    # find potential starts/ends of TADs
    det_scores = np.zeros([N,3],dtype=np.float64)
    det_scores = det_score(intgMAT2,np.arange(N),[W,W,W])
    S = np.sum(det_scores[:,0])/N
    det_scores[:,1] = det_scores[:,1]/S 
    det_scores[:,2] = det_scores[:,2]/S 

    start_peaks = np.asarray(detect_peaks(det_scores[:,2],  mpd=2),dtype=np.int64) #np.median(det_scores[:,1])
    end_peaks   = np.asarray(detect_peaks(det_scores[:,1],  mpd=2),dtype=np.int64) #np.median(det_scores[:,2])

    # we need to know the preceding end for each start
    prev_end = np.zeros(len(start_peaks),dtype=np.int64)
    k_old = 0
    for i in range(len(start_peaks)):
        k = k_old
        for j in range(k,len(end_peaks)):
            if end_peaks[k] <= start_peaks[i]:
                k = k + 1
            else:
                break
        if k==0:
            prev_end[i] = -1
            k_old = 0
        elif end_peaks[k-1] < start_peaks[0]:
            prev_end[i] = -1
            k_old = k-1
        else:
            prev_end[i] = k-1
            k_old = k - 1
    # we need to know the preceding start for each end
    prev_start = np.zeros(len(end_peaks),dtype=np.int64)
    k_old = 0
    for i in range(len(end_peaks)):
        k = k_old
        for j in range(k,len(start_peaks)):
            if start_peaks[k] <= end_peaks[i]:
                k = k + 1
            else:
                break
        if k==0:
            prev_start[i] = -1
            k_old = 0 
        else:
            prev_start[i] = k-1
            k_old = k - 1
    
    # pre-compute scores for start/end combinations
    M1 = len(start_peaks)
    M2 = len(end_peaks)

    rep = []
    indx = []
    for i in range(M1):
        ind = prev_end[i]
        indx.append(min(ind+1,len(end_peaks)))
        if ind == -1:
            rep.append(M2)
        elif ind == M2:
            rep.append(0)
        else:
            rep.append(M2 - ind - 1)
        
    long_i = np.repeat(start_peaks,rep)
    long_j = np.zeros(np.sum(rep),dtype=np.int64)   
    for i in range(M1):
        start = sum(rep[0:i])
        long_j[start:start+rep[i]] = end_peaks[indx[i]:len(end_peaks)] 

    long_l = long_j - long_i + 1
    temp_s = score(intgMAT2,long_i,long_l) #long_i is the center of a domain; long_l is the length of that domain

    # build a dense representation of the 'scores' but much smaller than N*N (for higher resolutions)
    i2I = np.zeros(N,dtype=np.int64)
    for i in range(len(start_peaks)):
        i2I[start_peaks[i]] = i
    j2J = np.zeros(N,dtype=np.int64)
    for j in range(len(end_peaks)):
        j2J[end_peaks[j]] = j

    scores = np.zeros([len(start_peaks),len(end_peaks)],dtype=np.float64)
    scores[i2I[long_i],j2J[long_j]] = temp_s

    # Dynamic Programming applies on potential start/end TAD boundaries we have found above
    T = np.zeros(len(end_peaks),dtype=np.float64)
    backT = np.zeros(len(end_peaks),dtype=np.int64) # shows the index of the start point (out of all start points) for each end point
    k_old = 0
    for j in range(len(end_peaks)):
        if prev_start[j] == -1:
            backT[j] = -1
            continue
        k = k_old
        kk = np.arange(int(prev_start[j])+1)  
        vals = parDP(scores,i2I,j2J,T,start_peaks,end_peaks,prev_end,j,kk)
        T[j] = np.max([np.max(vals),T[j-1]])
        backT[j] = np.argmax(vals)
        k_old = k

    # Extracting TADs based on backT
    counter = len(end_peaks) - 1
    tadCount = 0
    TADx1 = []
    TADx2 = []
    while True:# backT[counter] != 0 and backT[counter] != -1:
        
        if backT[counter] == -1:
            break
        elif backT[counter] == 0:
            TADx2.append(end_peaks[counter])
            TADx1.append(start_peaks[0])
            tadCount = tadCount + 1
            break
        elif counter > 0:
            TADx2.append(end_peaks[counter])
            TADx1.append(start_peaks[backT[counter]])
            tadCount = tadCount + 1
            counter = prev_end[backT[counter]]
            
            if prev_end[backT[counter]] == counter:
                counter = counter - 1
            if counter == -1:
                break
        elif counter == 0:
            TADx2.append(end_peaks[counter])
            TADx1.append(start_peaks[backT[counter]])
            tadCount = tadCount + 1
            break
    print('Number of TADs:',tadCount)
    TADx1.sort()
    TADx2.sort()
    TADx1 = (np.asarray(TADx1)).reshape([len(TADx1),1])
    TADx2 = (np.asarray(TADx2)).reshape([len(TADx2),1])
    TAD = np.concatenate((TADx1,TADx2),axis=1)
    np.savetxt(sys.argv[2],TAD,delimiter=' ',fmt='%d')

findTADs(intgMat)
