#!/usr/bin/env python

# Take 6-column BED and convert to 9 column GTF

import sys

def main(argv):

    bed = open(argv[1])
    gtf = open(argv[1].split(".bed")[0] + ".gtf", 'w')
    
    for line in bed:
        sline = line.split()
        chrom = sline[0]
        start = sline[1]
        end = sline[2]
        name = sline[3]
        score = sline[4]
        strand = sline[5]
        
        out_line = "{0}\tgene\tgene\t{1}\t{2}\t.\t{3}\t.\tgene_id {4}\n".format(chrom, start, end, strand, name)
        gtf.write(out_line)
        
    
if __name__ == '__main__':
    main(sys.argv)