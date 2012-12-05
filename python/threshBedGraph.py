#!/usr/bin/env python

import sys
import argparse

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', dest='bg', help='bedgraph')
    parser.add_argument('-t', dest='thresh', help="threshold")
    parser.add_argument('-o', dest='out_name', help='out track file')
    args = parser.parse_args()
    
    infile = open(args.bg)
    in_prefix = args.bg.split(".bed")[0]
    outfile = open(in_prefix + "_thresh" + args.thresh + ".bed", 'w')
    
    for line in infile:
        sline = line.split()
        if float(sline[3]) > float(args.thresh):
            outfile.write(line)
    
    infile.close()
    outfile.close()

if __name__ == '__main__':
    main(sys.argv)