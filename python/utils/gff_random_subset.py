#! /usr/bin/env python

# outputs random subset of records from a sorted GFF

import sys
import re
import random
import pdb

def main(argv):

    gff = open(argv[1])
    outfile = open(argv[2], 'w')
    fraction = float(argv[3])

    num_collected = 0
    num_total = 0
    gene = re.compile("gene")
    block=False
    #pdb.set_trace()
    for line in gff:
        sline = line.split()
        if gene.match(sline[2]):
            block=False
            num_total +=1
            if (random.random() <= fraction):
                num_collected += 1
                block=True
                outfile.write(line)
        elif block:
            outfile.write(line)

    print num_collected, num_total, float(num_collected)/num_total

if __name__ == "__main__":
    main(sys.argv)
