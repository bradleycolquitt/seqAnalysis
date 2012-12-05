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
    infile_prefix = args.input.split(".wig")[0]
    if not os.path.exists(infile_prefix): os.mkdir(infile_prefix)
    out_dir = infile_prefix
    print out_dir
    first = True
    outfile = ""
    for line in infile:
        if re.search("Step", line):
            if not first:
                outfile.close()
            sline = line.split()
            chr = sline[1].split("=")[1]
            outfile = open("/".join([out_dir, chr]) + ".wig", 'w')
            outfile.write(line)
            first = False
            continue
        outfile.write(line)
        
if __name__ == '__main__':
    main(sys.argv)