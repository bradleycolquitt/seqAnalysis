#!/usr/bin/env python

import sys
import os
import re
import argparse

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest="input")
    args = parser.parse_args()
    
    infile = open(args.input)
    
    outdir = args.input + "_chr"
    if not os.path.exists(outdir): os.mkdir(outdir)
    
    ind = range(1,20)
    ind.append('X')
    ind.append('Y')
    chrs = ["".join(['chr', str(chr)]) for chr in ind]
    
    file_dict = {}
    for chr in chrs:
        file_dict[chr] = open("/".join([outdir, chr]), 'w')
        
    for line in infile:
        sline = line.split()
        chr = sline[0]
        if chr in chrs:
            file_dict[chr].write(line)
        
if __name__ == '__main__':
    main(sys.argv)