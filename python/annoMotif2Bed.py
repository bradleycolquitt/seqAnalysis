#!/usr/bin/env python

import sys, os
import pdb

def main(argv):
    
    anno = open(argv[1])
    bed = open(argv[2], 'w')
    motif_length = int(argv[3])
    anno.readline()
    for line in anno:
        #pdb.set_trace()
        line = line.strip().split("\t")
        if len(line) < 21: continue
        chr = line[1]
        start = int(line[2])
        end = int(line[3])
        mid = (start + end) / 2
        offset = 0
        name = line[0]
        score = ""
        strand = "+"
        motif_info = line[20].split(")")
        try:
            if len(motif_info) > 2:
                count = 0
                for motif in motif_info:
                    if motif == '':
                        continue
                    count = count + 1
                    motif = motif.split("(")
                    #pdb.set_trace()
                    motif_0 = motif[0].split(",")
                    if len(motif_0) == 1:    
                        offset = int(motif_0[0])
                    elif len(motif_0) > 1:
                        offset = int(motif_0[1])
                    score = motif[1].split(",")[0]
                    strand = motif[1].split(",")[1]
                    start = mid + offset
                    end = start + motif_length
                    curr_name = name + "_" + str(count)
                    out = "\t".join([chr, str(start), str(end), curr_name, score, strand]) + "\n"
                    bed.write(out)
        except:
            pdb.set_trace()
if __name__ == '__main__':
    main(sys.argv)