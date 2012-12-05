#!/usr/bin/env python

import sys, re, pdb

# parse vista info file
# format:
#  element
#  hs position
#  hs id
#  mm position
#  mm id

def main(argv):
    
    infile = open(argv[1])
    outfile = open(argv[1].split(".")[0] + ".bed", 'w')
    
    name = ""
    chr = ""
    start = ""
    end = ""
    id = ""
    count = -1
    for line in infile:
        if re.search("mm", line) or re.search("hs", line):
            name = line.strip()
            count = 4
        if count == 1:  
            if re.search("chr", line):

                line1 = line.split(":")
                line2 = line1[1].strip().split("-")
                chr = line1[0]
                start = line2[0].replace(",", "")
                end = line2[1].replace(",", "")
        if count == 0:
            id = line.strip()
            name = "-".join([name, id])
        
            outline = "\t".join([chr, start, end, name, "0", "+"]) + "\n"
            outfile.write(outline)
            
        count -= 1    
    
if __name__ == '__main__':
    main(sys.argv)