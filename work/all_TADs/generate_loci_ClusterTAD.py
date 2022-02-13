import argparse
import numpy as np
import pandas as pd 

def main(file_in,file_out,chr,res):
    tad=np.loadtxt(file_in,skiprows=1,dtype=int)
    # print(tad)
    tad=tad[:,[1,3]]
    tad[:,0]=tad[:,0]-res
    # tad=tad.astype(int)
    # # print(tad)
    df_tad=pd.DataFrame(tad)
    df_tad[2]=chr
    df_tad[[2,0,1]].to_csv(file_out,sep='\t',header=False,index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', type=str, help='Input file armatus format')
    parser.add_argument('-o', type=str, help='Output file')
    parser.add_argument('-c', type=str, help='chromosome name')
    parser.add_argument('-r', type=int, help='Resolution')

    args = parser.parse_args()
    main(args.i, args.o, args.c, args.r)
    