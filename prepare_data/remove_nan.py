import numpy as  np  

def remove_nan(file_in,file_out):
    txt=np.loadtxt(file_in,dtype=np.float)
    txt[np.isnan(txt)]=0
    np.savetxt(file_out,txt,fmt='%g',delimiter='\t')

import sys, os

remove_nan(sys.argv[1],sys.argv[2])
if os.path.exists(sys.argv[1]):
    os.remove(sys.argv[1])
