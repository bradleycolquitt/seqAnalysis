#!/usr/bin/env python

import sys, os
import argparse
import pysam

def main(argv):
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', dest="bam")
    parser.add_argument('-t', dest="template")
    parser.add_argument('-o', dest="out")
    
    args = parser.parse_args()
    
    bam_in = pysam.Samfile(args.bam, 'rb')
    bam_template = pysam.Samfile(args.template, 'rb')
    bam_out = pysam.Samfile(args.out + ".bam", 'wb', template=bam_template)
    
    for read in bam_in.fetch():
        bam_out.write(read)
        
    bam_in.close()
    bam_template.close()
    bam_out.close()
    
    bam_sorted = bam_out + "_sort.bam"
    pysam.sort(bam_out + ".bam", bam_sorted)
    pysam.index(bam_sorted)
    
if __name__ == '__main__':
    main(sys.argv)
