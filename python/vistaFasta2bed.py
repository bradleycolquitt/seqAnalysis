#!/usr/bin/env python

import sys, re

def main(argv):
    
    infile = open(argv[1])
    outfile = open(argv[1].split(".fa")[0] + ".bed", 'w')
    
    expr = ""
    for line in infile:
        if re.search(">", line):
            line = line.split("|")
            pos = line[1].strip()
            pos1 = pos.split(":")
            chr = pos1[0]
            (start, end) = pos1[1].split("-")
            name = "-".join([line[0].split(">")[1], "_".join(line[2].strip().split())])
                
            type = "1"
            if line[3].strip() == "negative": type = "0"
            
            if len(line) == 5:
                expr = line[4].split("[")[0]
            else:
                expr = ""
            outline = "\t".join([chr, start, end, name, type, "+", expr]) + "\n"
            outfile.write(outline)
    
    
if __name__ == '__main__':
    main(sys.argv)