import argparse
import numpy as np

def main(file_in,file_out,chr):
    ff = open(file_out, 'w')
    with open(file_in, 'r') as f:
        line = f.readlines()
        count = len(line)
        for cnt in range(count):
            thisline = line[cnt]
            index1 = thisline.find('\t')
            if cnt == (count - 1):
                index2 = thisline.find('\n', index1 + 1)
                #print(index2)
                if index2 != -1:
                    ff.write(str(chr) + '\t' + thisline[index1 + 1:])
            else:
                ff.write(str(chr) + '\t' + thisline[index1 + 1:])
    ff.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', type=str, help='Input file DI format')
    parser.add_argument('-o', type=str, help='Output file')
    parser.add_argument('-c', type=str, help='Chromosome name')

    args = parser.parse_args()
    main(args.i, args.o, args.c)
