import numpy as np 
import os
import sys

def convert(file_in,file_out):
    mat=np.loadtxt(file_in)
    # print(mat)
    with open(file_out,'w') as f :
        for i in range(mat.shape[0]):
            for j in range(i, mat.shape[1]):
                if mat[i,j]!=0:
                    # print(mat[i,j])
                    f.write('%g\t%g\t%g\n'%((i+1),(j+1),mat[i,j]))

convert(sys.argv[1],sys.argv[2])

