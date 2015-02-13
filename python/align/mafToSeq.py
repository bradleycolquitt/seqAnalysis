#! /usr/bin/env python

"""Extract query sequence from a MAF file. Expects MAF to be sorted by target"""

import sys
import os.path
import re
import fasta_utils as futil
import pdb
def main(argv):

    infile = open(argv[1])
    outfile = open(os.path.splitext(argv[1])[0] + ".fa", 'w', -1)

    pos = 0
    curr_target_name = ""
    start = 0
    a_pattern = re.compile("^a")
    s_pattern = re.compile("^s")
    for line in infile:
        if s_pattern.search(line):
            sline = line.split()
            if pos == 0:
                # target
                pos = 1
                target = sline[1]
                if curr_target_name != target:
                    if curr_target_name != "":
                        outfile.write("\n")
                        start = 0
                    outfile.write(">{0}\n".format(target))
                    curr_target_name = target
            elif pos == 1:
                # query
                pos = 0
                seq = sline[6].translate(None, "-")
                if sline[4] == "-":
                    seq = futil.reverse_complement(seq)
                start = futil.write_segment(seq, start, outfile)
    outfile.write("\n")

if __name__ == "__main__":
    main(sys.argv)
