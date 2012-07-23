#!/usr/bin/env python

import sys
import argparse
import itertools
import re
import numpy as np

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', dest='one_name', help='wig file')
    #parser.add_argument('-b', dest='two_name', help='wig file')
    parser.add_argument('-t', dest='thresh', help="threshold value", type=float)
    parser.add_argument('-o', dest='out_name', help='wig file')    
    args = parser.parse_args()
    
    one = open(args.one_name)
#    two = open(args.two_name)
    out = open(args.out_name, 'w')
    
    for line in one:
        if re.search("chr", line):
            #assert line1 == line2, "line1 is %s, line2 is %s" % (line1,line2)
            out.write(line)
            continue
        line = float(line.strip())
        #line2 = float(line2.strip())
        out_val = 0
        if line > args.thresh: out_val = 1
        out.write(str(out_val) + "\n")        
    out.close()
    
if __name__ == '__main__':
    main(sys.argv)
