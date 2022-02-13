import argparse
import numpy as np
import pandas as pd 

def main(file_in,file_out,res):
    tad=pd.read_csv(file_in,header=None,usecols=(1,2),sep='\t',dtype=int)
    tad.loc[:,2]=tad.loc[:,2]+1
    tad.loc[:,1]=tad.loc[:,1]+res
    tad=tad/res
    tad=tad.astype(int)
    # print(tad)
    tad.to_csv(file_out,sep='\t',header=False,index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', type=str, help='Input file armatus format')
    parser.add_argument('-o', type=str, help='Output file')
    parser.add_argument('-r', type=int, help='Resolution')

    args = parser.parse_args()
    main(args.i, args.o, args.r)
    
