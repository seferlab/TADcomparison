import argparse
import numpy as np
import pandas as pd 

def main(file_in,file_out):
    tad=np.loadtxt(file_in,dtype=int,usecols=(0,1),skiprows=1)
    tad=tad+1
    np.savetxt(file_out,tad,fmt='%g',delimiter='\t')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', type=str, help='Input file armatus format')
    parser.add_argument('-o', type=str, help='Output file')

    args = parser.parse_args()
    main(args.i, args.o)
    