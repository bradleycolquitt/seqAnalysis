#! /usr/bin/env python

import os, sys
from string import *
from itertools import izip

def main(argv):
    a = open(argv[1], 'r')
    b = open(argv[2], 'r')
    a_name = argv[1].split(".bed")[0]
    b_name = argv[2].split(".bed")[0]
    outfile_name = a_name + "_sub_" + b_name
    outfile = open(outfile_name, 'w')
    with a as f1:
        with b as f2:
            for (l1, l2) in izip(f1, f2):
                ls1 = l1.strip().split()
                ls2 = l2.strip().split()
                val = atoi(ls1[4]) - atoi(ls2[4])
                if val < 0:
                    val = 0
                out = "\t".join(ls1[0:4]) + "\t"+ str(val) + "\n"
                outfile.write(out)
    outfile.close()
if __name__ == "__main__":
    main(sys.argv)
