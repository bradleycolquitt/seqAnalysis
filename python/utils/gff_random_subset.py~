#! /usr/bin/env python

# splits sorted GFF into individual genes

import sys
import re

def main(argv):

    gff = open(argv[1])
    outdir = argv[2]

    curr_file = 0
    gene = re.compile("gene")
    for line in gff:
        sline = line.split()
        if gene.match(sline[2]):
            outname = sline[8].split(";")[0].split("=")[1]
            if curr_file > 0:
                curr_file.close()
            curr_file = open(outdir + outname + ".gff", 'w')
        curr_file.write(line)


if __name__ == "__main__":
    main(sys.argv)
