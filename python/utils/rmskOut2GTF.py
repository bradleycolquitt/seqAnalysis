#! /usr/bin/env python

import sys
import re
import pdb
def main(argv):

    infile = open(argv[1])
    outfile = open(argv[1] + ".gtf", 'w')
    #pdb.set_trace()
    num = re.compile("[0-9]+")
    for line in infile:
        sline = line.split()
        if len(sline) > 1:
            if num.match(sline[0]):
                out = '{0}\trmsk\t{1}\t{2}\t{3}\t.\t+\t.\tgene_id "{4}-{5}-{6}-{7}"\n'.format(
                    sline[4],
                    sline[10],
                    sline[5],
                    sline[6],
                    sline[14],
                    sline[11],
                    sline[12],
                    sline[13]
                )
                outfile.write(out)

if __name__ == "__main__":
    main(sys.argv)
