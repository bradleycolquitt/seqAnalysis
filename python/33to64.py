#!/usr/bin/env python

import sys
import re

def main(argv):
    
    infile = open(argv[1])
    outfile = open(argv[1] + "_64", 'w')
    qual = False
    for line in infile:
        if qual:
            line = line.strip()
            nline = []
            for char in list(line):
                nline.append(chr(ord(char) + 31))
            nline = "".join(nline) + "\n"
            outfile.write(nline)
            qual = False
            continue    
        pat = re.compile("^\+$")
        if re.search(pat, line):
            qual = True
        outfile.write(line)
    infile.close()
    outfile.close()
    
if __name__ == '__main__':
    main(sys.argv)