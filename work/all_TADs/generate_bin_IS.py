import argparse
import numpy as np
import pandas as pd 

def main(file_in,file_out,chr_length):
    tad=np.loadtxt(file_in,usecols=(5),skiprows=1)
    #print(tad)
    tad=tad.astype(int)
    #print(tad)
    final_tad=np.empty(shape=[0, 2],dtype=int)
    stt=1
    len_boundary=len(tad)
    for i in range(len_boundary):
        if i!=len_boundary-1:
            final_tad=np.append(final_tad, [[stt,tad[i]]], axis=0)
            stt=tad[i]+1
        else:
            final_tad=np.append(final_tad, [[stt,tad[i]]], axis=0)
            final_tad=np.append(final_tad, [[tad[i]+1,chr_length]], axis=0)
    # print(final_tad)
    np.savetxt(file_out,final_tad,fmt='%g',delimiter='\t')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', type=str, help='Input file armatus format')
    parser.add_argument('-o', type=str, help='Output file')
    parser.add_argument('-l', type=int, help='chromosome length')

    args = parser.parse_args()
    main(args.i, args.o, args.l)
