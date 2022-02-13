import sys, getopt
import numpy as np

def main(argv):
   inputfile = ''
   outputfile = ''
   chrname=''
   binsize=0
   try:
      opts, args = getopt.getopt(argv,"hi:o:c:b:",["ifile=","ofile=","chrom_name","binsize"])
   except getopt.GetoptError:
      print('test.py -i <inputfile> -o <outputfile> -c <chromosome name> -b <binsize>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('test.py -i <inputfile> -o <outputfile> -c <chromosome name> -b <binsize>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
      elif opt in ("-c", "--chrom_name"):
          chrname = arg
      elif opt in ("-b", "--binsize"):
          binsize = int(arg)
   # print('input:', inputfile)
   # print('output:', outputfile)
   ff = open(outputfile, 'w')
   mat = np.loadtxt(inputfile).astype(np.float)
   n_rows = mat.shape[0]
   for i in range(0, n_rows):
       for j in range(i, n_rows):
           if mat[i][j] != 0:
               line_new='0'+'\t'+str(chrname)+'\t'+str(i*binsize)+'\t'+'0'+'\t'+'0'+'\t'+str(chrname)+'\t'+str(j*binsize)+'\t'+'1'+'\t'+str(mat[i][j])+'\n'
               ff.write(line_new)
   ff.close()

if __name__ == "__main__":
   main(sys.argv[1:])
