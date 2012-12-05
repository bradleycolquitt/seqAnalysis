#! /usr/bin/env python

import sys
from string import *

def main(argv):
    infile_name = argv[1].split(".bed")[0]
    infile = open(argv[1])
    thresh = atoi(argv[2])
    outfile = open(infile_name + "_thresh" + argv[2] + ".bed", 'w')    
    for line in infile:
        sline = line.strip().split()
        line_chr = sline[0]
        line_start = atoi(sline[1])
        line_end = atoi(sline[2])
        span = line_end - line_start
        if span >= thresh:
            out = line_chr + "\t" + str(line_start) + "\t" + str(line_end) + "\n"
            outfile.write(out)
    outfile.close()

if __name__ == "__main__":
    main(sys.argv)
