#!/usr/bin/env python

import sys

def main(argv):
    
    infile = open(argv[1])
    outfile = open(argv[1] + ".bed", 'w')
    line = infile.readline()
    
    for line in infile:
        line = line.split()
        start = line[18].split(")")[0]
        end = line[19].split("(")[1].split(",")[0]
        out = "\t".join([line[15], start, end, line[0], line[8], line[16]]) + "\n"
        outfile.write(out)
    
if __name__ == '__main__':
    main(sys.argv)