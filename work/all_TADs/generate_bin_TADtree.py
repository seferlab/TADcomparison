import argparse
import numpy as np
import pandas as pd 

def main(file_in,file_out):
    tad=pd.read_csv(file_in,skiprows=1,header=None,sep='\t',usecols=(1,2),dtype=int)
    tad=tad+1
    # print(tad)
    tad.to_csv(file_out,sep='\t',header=False,index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', type=str, help='Input file tadtree format')
    parser.add_argument('-o', type=str, help='Output file')

    args = parser.parse_args()
    main(args.i, args.o)
