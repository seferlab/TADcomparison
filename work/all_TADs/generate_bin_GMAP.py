import argparse
import numpy as np
import pandas as pd 

def main(file_in,file_out,res):
    res=int(res)
    mid_res=int(res/2)
    tad=pd.read_csv(file_in,skiprows=0,header=None,usecols=(0,1),sep='\t',dtype=int)
    tad=tad+mid_res
    tad.iat[0,0]=tad.iat[0,0]-res
    tad=(tad/res).astype(np.int)
    tad[[0]]=tad[[0]]+1
    tad[[0,1]].to_csv(file_out,sep='\t',header=False,index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', type=str, help='Input file tadtree format')
    parser.add_argument('-o', type=str, help='Output file')
    parser.add_argument('-r', type=int, help='Resolution')

    args = parser.parse_args()
    main(args.i, args.o, args.r)
