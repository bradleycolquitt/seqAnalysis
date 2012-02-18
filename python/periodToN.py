#!/usr/bin/env python

import sys
import string

def main(argv):
    #print argv
    infile = open(argv[1])
    outfile = open(argv[1] + "_Nconvert", 'w')
    #qual = False
    seq = False
    for line in infile:
        if seq:
            line = line.replace(".", "N")
            outfile.write(line)
            seq = False
            continue
        #sline = line.split()
        #ssline = list(sline[0])    
        if "@QSEQ" in line:
            seq = True
        #if ssline[0] == "+":
        #    qual = True
            
        outfile.write(line)
    infile.close()
    outfile.close()
    
if __name__ == '__main__':
    main(sys.argv)