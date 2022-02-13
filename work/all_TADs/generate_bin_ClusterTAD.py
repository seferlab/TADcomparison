import argparse
import numpy as np
import pandas as pd 

def main(file_in,file_out):
    tad=np.loadtxt(file_in,skiprows=1,dtype=int)
    # print(tad)
    tad=tad[:,[0,2]]
    # tad=tad.astype(int)
    # # print(tad)
    np.savetxt(file_out,tad,delimiter='\t',fmt='%d')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', type=str, help='Input file armatus format')
    parser.add_argument('-o', type=str, help='Output file')

    args = parser.parse_args()
    main(args.i, args.o)
    