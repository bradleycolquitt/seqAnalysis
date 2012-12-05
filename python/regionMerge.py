#! /usr/bin/env python

import sys
from string import *

def main(argv):
    infile = open(argv[1])
    gap = atoi(argv[2])
    outfile = open(argv[1] + "_merge" + argv[2], 'w')

    chr_store = ""
    start_store = 0
    end_store = 0
    
    
    for line in infile:
        sline = line.strip().split()
        line_chr = sline[0]
        line_start = atoi(sline[1])
        line_end = atoi(sline[2])
        if end_store == 0:
            chr_store = line_chr
            start_store = line_start
            end_store = line_end
            continue
        gap_pos = end_store + gap
        if line_start > gap_pos and chr_store == line_chr:
            out = sline[0] + "\t" + str(start_store) + "\t" + str(end_store) + "\n"
            outfile.write(out)
            start_store = line_start
            end_store = line_end
        elif chr_store != line_chr:
            end_store = 0
        else:
            end_store = atoi(sline[2])
    
    outfile.close()

if __name__ == "__main__":
    main(sys.argv)
