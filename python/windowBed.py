#!/usr/bin/env python

import sys
import os
import argparse
import pysam
from multiprocessing import Queue, Process

chrs = []
for i in [range(19) + 1, 'X', 'Y', 'M']:
    chrs.append("".join(['chr', i]))
    
class windower:
    def __init__(self, bamfile, outpath, window):
        self.bamfile = bamfile
        self.outpath = outpath
        self.window = window
        self.chr_lengths = bamfile.lengths
        #for line in genome:
	#    line = line.split()
	#    if len(line)>0:
        #        chr_lengths[line[0]] = line[1]
        self.chrs_queue = Queue()
        for chr in bamfile.references:
            q.put(chr)
    
    def window(self):
        ## Setup interface to parallel proc
        ## Iterate through chrs_queue
        p = Process(target=window_core, args=(q,out_path))
        p.start()
        p.join()
        
    def window_core(chr, out_path):
        bedfile = open(out_path + "/" + chr, 'w')
        chr_length = self.chr_lengths(chr)
        curr_pos = 0
        next_pos = curr_pos + self.window
        count = 0
        while curr_pos <= chr_length:
            for read in self.bamfile.fetch(chr):
                read_mid = read.pos + read.isize / 2
                if read_mid > curr_pos and read_mid < next_pos:
                    count = count + 1
                elif read_mid > next_pos:
                    out = "/t".join([chr, str(curr_pos), str(count)])
                    bedfile.write(out)
        bedfile.close()

def main(argv):
    ## Take BAM of reads
    ## Count number of reads within specified window size
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', dest='bam')
    parser.add_argument('-w', dest='window')
    args = parser.parse_args()
    genome = open("/seq/lib/mouse.mm9.genome", 'r')
                 
    bamfile = pysam.Samfile(args.bam, 'rb')
    outpath = "/".join([bedpath, args.bam.split(".bam")[0]])
    w = windower(bamfile, outpath, args.window)
if __name__ == "__main__":
    main(sys.argv)