import argparse
import numpy as np

def main(file_in,file_out):
    tad=np.loadtxt(file_in,dtype=int)
    # print(tad)
    final_tad=np.empty(shape=[0, 2],dtype=int)
    stt=1
    for i in range(len(tad)):
        final_tad=np.append(final_tad, [[stt,tad[i]]], axis=0)
        stt=tad[i]+1
    np.savetxt(file_out,final_tad,fmt='%g',delimiter='\t')




if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', type=str, help='Input file hicseg format')
    parser.add_argument('-o', type=str, help='Output file')

    args = parser.parse_args()
    main(args.i, args.o)
