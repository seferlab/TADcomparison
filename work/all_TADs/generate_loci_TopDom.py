import argparse
import numpy as np
import pandas as pd

def main(file_in,file_out):
    tad=pd.read_csv(file_in,header=None,sep='\t')
    tad=tad[tad[3]=='domain']
    tad=tad[[0,1,2]]
    # tad=tad.astype(int)
    # print(tad)
    tad.to_csv(file_out,sep='\t',header=False,index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', type=str, help='Input file topdom format')
    parser.add_argument('-o', type=str, help='Output file')

    args = parser.parse_args()
    main(args.i, args.o)
