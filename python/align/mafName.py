#! /usr/bin/env python

"""
Adds species names to a MAF file

Usage: mafName.py file.maf src1 src2

Write new file "_named.maf"
"""

import sys
import os
import re

def main(argv):
    ifname = argv[1]
    name1 = argv[2]
    name2 = argv[3]

    ofname = os.path.splitext(ifname)[0] + "_named.maf"

    ifile = open(ifname)
    ofile = open(ofname, "w")

    pos = 0
    a_pattern = re.compile("^a")
    s_pattern = re.compile("^s")
    for line in ifile:
        if s_pattern.search(line):
            sline = line.split()
            try:
                if pos == 0:
                    # target
                    sline[1] = ".".join([name1, sline[1]])
                    ofile.write("\t".join(sline) + "\n")
                    pos = 1
                elif pos == 1:
                    # query
                    sline[1] = ".".join([name2, sline[1]])
                    ofile.write("\t".join(sline) + "\n")
                    pos = 0
            except:
                print "Malformed line: " line
        else:
            ofile.write(line)

if __name__ == "__main__":
    main(sys.argv)
