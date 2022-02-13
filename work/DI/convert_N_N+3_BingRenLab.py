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
   cnt = 0
   with open(inputfile, 'r') as f:
       line = f.readlines()
       for line_list in line:
           if (line_list != ""):
               binstart = cnt * binsize
               binend = (cnt+1) * binsize
               line_new = chrname+ '\t' + str(binstart) + '\t' + str(binend) + line_list
               ff.write(line_new)
               cnt = cnt + 1
   ff.close()

if __name__ == "__main__":
   main(sys.argv[1:])