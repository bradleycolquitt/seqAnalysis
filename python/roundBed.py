#!/usr/bin/env python

import sys
from string import *
from math import *

def rounder(x,y):
    return str(int(round(atoi(x),y)))

def main(argv):
    
    round_factor = -1 * int(log10(atoi(argv[1])))
    input_name = argv[2]
    infile = open(input_name, 'r')
    outfile = open(input_name + "_round" + argv[1], 'w')
    for line in infile:
        line = line.strip()
        sline = line.split()
        out = sline[0] + "\t" + rounder(sline[1], round_factor) + "\t" + rounder(sline[2], round_factor) + "\t" + "\t".join(sline[3:]) + "\n"
        outfile.write(out)
    infile.close()
    outfile.close()

if __name__ == "__main__":
    main(sys.argv)

