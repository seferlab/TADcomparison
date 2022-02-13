import argparse
import numpy as np

def main(file_in,file_out):
    ff = open(file_out, 'w')
    with open(file_in, 'r') as f:
        end_mat = np.loadtxt(file_in, usecols=(2), dtype=int)
        line = f.readlines()
        count = len(line)
        for cnt in range(count - 1, -1, -1):
            thisline = line[cnt]
            if (thisline != ""):
                index1 = thisline.find('\t')
                index2 = thisline.find('\t', index1 + 1)
                ff.write('chr' + thisline[:index2] + '\t' + str(end_mat[cnt] + 1) + '\n')
    ff.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', type=str, help='Input file armatus format')
    parser.add_argument('-o', type=str, help='Output file')

    args = parser.parse_args()
    main(args.i, args.o)
