#! /usr/bin/env python

import sys
import os
import pysam
import pdb

cpdef thresh(ifname, threshin):
    cdef str ofname = os.path.splitext(os.path.basename(ifname))[0] + "_C" + threshin + ".bam"
    bamfile = pysam.AlignmentFile(ifname)
    outfile = pysam.AlignmentFile(ofname, 'wb', bamfile)

    cdef int thresh = int(threshin)
    cdef int count = 0
    cdef int pos = -1
    cdef int curr_pos = 0
    cdef object read
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

cpdef main(argv):

    cdef str ifname = argv[1]
    thresh(ifname, argv[2])



if __name__ == "__main__":
   main(sys.argv)
