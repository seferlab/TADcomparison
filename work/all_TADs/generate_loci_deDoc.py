import argparse
import numpy as np

def main(file_in,file_out,chr,res):
    
    ff = open(file_out, 'w')
    numfile=open(file_in+'(M)')
    numlist=[]
    for line in numfile:
        line=line.strip()
        numlist=line.split(' ')
        start_bin=int(numlist[0])-1
        end_bin=int(numlist[-1])
        # print(numlist[0])
        # print(numlist[-1])
        ff.write(str(chr) + '\t' + str(start_bin * res) + '\t' + str( end_bin * res) + '\n')
    numfile=open(file_in+'(E)')
    numlist=[]
    for line in numfile:
        line=line.strip()
        numlist=line.split(' ')
        start_bin=int(numlist[0])-1
        end_bin=int(numlist[-1])
        # print(numlist[0])
        # print(numlist[-1])
        ff.write(str(chr) + '\t' + str(start_bin * res) + '\t' + str( end_bin * res) + '\n')
    ff.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', type=str, help='Input file tadtree format')
    parser.add_argument('-o', type=str, help='Output file')
    parser.add_argument('-c', type=str, help='Chromosome name')
    parser.add_argument('-r', type=int, help='Resolution')

    args = parser.parse_args()
    main(args.i, args.o, args.c, args.r)

