import sys, getopt
from itertools import islice

def main(argv):
    inputfile = ''
    outputfile = ''
    try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
      print('test.py -i <inputfile> -o <outputfile>')
      sys.exit(2)
    for opt, arg in opts:
      if opt == '-h':
         print('test.py -i <inputfile> -o <outputfile>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
    # print('input:', inputfile)
    # print('output:', outputfile)
    ff = open(outputfile, 'w')
    with open(inputfile, 'r') as f:
        for line in islice(f, 1, None):
            pos = line.find('\t')
            print(pos)
            ff.write(line[pos + 1:])
    ff.close()

if __name__ == "__main__":
   main(sys.argv[1:])