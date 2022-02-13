import sys
import numpy as np

def main(file_out1,file_out2,bin_size,mat_dim,chr_name):
    bin_size=int(bin_size)
    mat_dim=int(mat_dim)
    if chr_name=='chrX':
        chr_idx=23
    else:
        chr_idx=chr_name[3:]
        chr_idx=int(chr_idx)
    filename = file_out1
    with open(filename, 'w') as f:
        for i in range(1,(chr_idx+1)):
            if i==chr_idx:
                f.write("%d\t%s\t1\t%d\n"%(chr_idx,chr_name,mat_dim) )
            else:
                f.write("%d\tchr%d\t0\t0\n"%(i,i) )
    filename = file_out2
    with open(filename, 'w') as f:
        for i in range(int(mat_dim)):
            f.write("%d\t%d\t%d\n"%((chr_idx-1),(i*bin_size+1),(i+1)*bin_size) )

if __name__ == "__main__":
   main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])

