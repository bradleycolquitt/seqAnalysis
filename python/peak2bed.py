#!/usr/bin/env python

import sys

def main(argv):

    infile = open(argv[1])
    outfile = open(argv[1].split(".")[0] + ".bed", 'w')
    
    for line in infile:
        line = line.split()
        
        outline = "\t".join(line[0:3] + line[6:8] + ["+"]) + "\n"
        outfile.write(outline)
    
if __name__ == '__main__':
    main(sys.argv)