#!/usr/bin/env python

import sys
import os
import re
import pysam

def main(argv):
    read_ids = []
    mate1 = pysam.Samfile(argv[1], 'r')
    case = 1
    mate2 = ""
    mate2_name = ""
    out = ""
    if re.search("fastq", argv[2]):
        case = 1
        mate2 = open(argv[2])
        mate2_name = os.path.basename(argv[2])
        out = open("_".join([argv[1], mate2_name]), 'w')
    else:
        case = 2
        mate2 = pysam.Samfile(argv[2], 'r')
        mate2_name = os.path.basename(argv[2])
        out = pysam.Samfile("_".join([argv[1], mate2_name]), 'w', template=mate2)
    for read in mate1:
        read_ids.append(read.qname)
    ind = 0
    if case == 1:
        while True:
            read = mate2.readline()
            if read == "": break
            if re.search("@", read):
                id = read.split()[0]
                id = id.split('@')[1]
                if id in read_ids:
                    out.write(read)
                    while ind < 4:
                        out.write(mate2.readline())
                        ind = ind + 1
                    ind = 0
    else:
        for read in mate2:
            if read.qname in read_ids:
                out.write(read)
            
        
    
if __name__ == '__main__':
    main(sys.argv)