import argparse
import numpy as np

def main(file_in,file_out,chr,chr_len,res):
    ff = open(file_out, 'w')
    mat = np.loadtxt(file_in, skiprows=1, usecols=5).astype(int)
    bound_list = mat.tolist()
    cnts = len(bound_list)
    for i in range(cnts):
        if i == 0:
            ff.write(str(chr) + '\t0\t' + str(bound_list[i] * res) + '\n')
        elif i == cnts - 1:
            ff.write(str(chr) + '\t' + str(bound_list[i - 1] * res) + '\t' + str(chr_len * res) + '\n')
        else:
            ff.write(str(chr) + '\t' + str(bound_list[i - 1] * res) + '\t' + str(bound_list[i] * res) + '\n')
    ff.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', type=str, help='Input file insulation format')
    parser.add_argument('-o', type=str, help='Output file')
    parser.add_argument('-c', type=str, help='Chromosome name')
    parser.add_argument('-l', type=int, help='Chromosome length')
    parser.add_argument('-r', type=int, help='Resolution')

    args = parser.parse_args()
    main(args.i, args.o, args.c, args.l, args.r)
