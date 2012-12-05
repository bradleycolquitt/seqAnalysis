#!/usr/bin/env python

import sys, re

def main(argv):
    
    infile = open(argv[1])
    outfile = open(argv[1] + "_rm", 'w')
    
    ind = range(1,20)
    ind.append('X')
    ind.append('Y')
    chrs = ["".join(['chr', str(chr)]) for chr in ind]
    print chrs
    cont = False
    for line in infile:
        if re.search("fixedStep", line):
            #print line
            cont = False
            sline = line.strip().split()
            chr = sline[1].split("=")[1]
            if chr in chrs:
                print chr
                outfile.write(line)
                cont = True
        elif cont:
            outfile.write(line)
    
    infile.close()
    outfile.close()
    
            
if __name__ == '__main__':
    main(sys.argv)
