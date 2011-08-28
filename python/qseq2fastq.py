#!/usr/bin/env python

import sys

def main(argv):

    infile = open(argv[1], 'r')
    outfile = open(argv[2], 'w')

    for line in infile:
        line = line.strip()
        sline = line.split()

        out = "@" + ":".join(sline[0:5]) + "#" + sline[6] + "/" + sline[7] + "\n" \
              + sline[8] + "\n" \
              + "+" + sline[0] + ":" + ":".join(sline[2:5]) + "#" + sline[6] + "/" + sline[7] + "\n" \
              + sline[9] + "\n"
        outfile.write(out)
    
    infile.close()
    outfile.close()

if __name__ == "__main__":
    main(sys.argv)
