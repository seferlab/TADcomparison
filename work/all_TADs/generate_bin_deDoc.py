import argparse
import numpy as np

def main(file_in,file_out):
    ff = open(file_out, 'w')
    numfile=open(file_in+'(M)')
    numlist=[]
    for line in numfile:
        line=line.strip()
        numlist=line.split(' ')
        start_bin=int(numlist[0])
        end_bin=int(numlist[-1])
        # print(numlist[0])
        # print(numlist[-1])
        ff.write(str(start_bin) + '\t' + str( end_bin) + '\n')
    numfile=open(file_in+'(E)')
    numlist=[]
    for line in numfile:
        line=line.strip()
        numlist=line.split(' ')
        start_bin=int(numlist[0])
        end_bin=int(numlist[-1])
        # print(numlist[0])
        # print(numlist[-1])
        ff.write(str(start_bin) + '\t' + str( end_bin) + '\n')
    ff.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', type=str, help='Input file tadtree format')
    parser.add_argument('-o', type=str, help='Output file')

    args = parser.parse_args()
    main(args.i, args.o)

