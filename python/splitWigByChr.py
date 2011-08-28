#!/usr/bin/env python

import sys
import re
import argparse

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest="input")
    args = parser.parse_args()
    
    infile = open(args.input)
    first = True
    outfile = ""
    for line in infile:
        if re.search("Step", line):
            if not first:
                outfile.close()
            sline = line.split()
            chr = sline[2].split("=")[1]
            outfile = open("/".join([out_dir, chr]), 'w')
            outfile.write(line)
            first = False
        outfile.write(line)
        
if __name__ == '__main__':
    main(sys.argv)