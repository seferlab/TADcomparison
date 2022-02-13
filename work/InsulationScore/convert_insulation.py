import sys, getopt

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
    with open(inputfile, 'r') as f:
        line = f.readlines()
        line_header=''
        for counter in range(0,len(line)):
            binstart=counter*binsize+1
            binend = (counter + 1) * binsize + 1
            line_header = line_header+'\t'+'bin'+ str(counter) + '|ce40|' + str(chrname) + ':' + str(binstart) + '-' + str(binend)
        line_header=line_header+'\n'
        ff.write(line_header)
        cnt = 0
        for line_list in line:
            if(line_list !=""):
                binstart = cnt * binsize+1
                binend = (cnt + 1) * binsize+1
                line_new ='bin'+ str(cnt)+ '|ce40|' + str(chrname) + ':' + str(binstart) + '-' + str(binend) +'\t'+ line_list
                ff.write(line_new)
                cnt = cnt + 1
    ff.close()


if __name__ == "__main__":
    main(sys.argv[1:])