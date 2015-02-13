#! /usr/bin/env python

import sys
import os.path
import re
import pdb

def main(argv):

    infile = open(argv[1])
    out_fname = os.path.splitext(argv[1])[0]

    out_files = {}
    #out_target = open( + ".bed", 'w')
    #out_query = open(os.path.splitext(argv[1])[0] + ".bed", 'w')

    score = ""
    ind = 0
    a_pattern = re.compile("^a")
    s_pattern = re.compile("^s")
    for line in infile:
        if a_pattern.search(line):
            #pdb.set_trace()
            score = line.split()[1].split("=")[1].split(".")[0]
            ind += 1
        elif s_pattern.search(line):
            #pdb.set_trace()
            sline = line.split()
            sline1s = sline[1].split(".")

            outline = "\t".join([sline1s[1],
                                 str(int(sline[2]) ),
                                 str(int(sline[2]) + int(sline[3])),
                                 "maf"+str(ind),
                                 "0",
                                 sline[4]]) + "\n"
            if sline1s[0] in out_files:
                ret = out_files[sline1s[0]].write(outline)
            else:
                out_files[sline1s[0]] = open(out_fname + "_" + sline1s[0] + ".bed", 'w')

if __name__ == "__main__":
    main(sys.argv)
