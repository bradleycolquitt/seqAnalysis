#! /usr/bin/env python

import sys
import os
import pysam
import pdb

def thresh(ifname, threshin):
    ofname = os.path.splitext(os.path.basename(ifname))[0] + "_C" + threshin + ".bam"
    bamfile = pysam.AlignmentFile(ifname)
    outfile = pysam.AlignmentFile(ofname, 'wb', bamfile)

    thresh = int(threshin)
    count = 0
    pos = -1
    curr_pos = 0
    #pdb.set_trace()
    for read in bamfile.fetch():
        if count > thresh: continue
        outfile.write(read)

        # if read.is_reverse:
        #    curr_pos = read.alen
        # else:
        #    curr_pos = read.pos

        curr_pos = read.pos
        if curr_pos == pos:
           count += 1
        else:
           pos = curr_pos
           count = 0

def main(argv):

    ifname = argv[1]
    thresh(ifname, argv[2])



if __name__ == "__main__":
   main(sys.argv)
