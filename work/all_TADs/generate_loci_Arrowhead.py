import argparse
import numpy as np

def main(file_in,file_out):
    ff = open(file_out, 'w')
    with open(file_in, 'r') as f:
        line = f.readlines()
        count = len(line)
        for cnt in range(2, count):
            thisline = line[cnt]
            if (thisline != ""):
                index1 = thisline.find('\t')
                index2 = thisline.find('\t', index1 + 1)
                index3 = thisline.find('\t', index2 + 1)
                ff.write('chr' + thisline[:index3] + '\n')
    ff.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', type=str, help='Input file arrowhead format')
    parser.add_argument('-o', type=str, help='Output file')

    args = parser.parse_args()
    main(args.i, args.o)
