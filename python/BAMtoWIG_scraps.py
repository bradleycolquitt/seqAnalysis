#!/usr/bin/env python

import sys
import os
import argparse
import pysam
import numpy as np
from multiprocessing import Queue, Process, Pool
from string import *

wigpath = "/media/storage2/data/wig"

class windower:
    def __init__(self, bamname, outpath, window_size):
        self.bamname = bamname
        self.bamfile = pysam.Samfile(bamname, 'rb')
        self.outpath = outpath
        self.window_size = atoi(window_size)
        self.nreads = atoi(pysam.flagstat(bamname)[0].split()[0])
        self.chr_lengths = self.bamfile.lengths
        self.chrs_queue = []
        for index in range(self.bamfile.nreferences):
            self.chrs_queue.append((self.bamfile.references[index], self.bamfile.lengths[index]))
    
    def window(self):
        #print self.chrs_queue
        #sys.exit()
        p = Pool(processes=2)
        args = []
        for chr in self.chrs_queue:
            #window_core(self.bamfile, self.window_size, chr, self.nreads, self.outpath)
            args.append((self.bamname, self.window_size, chr, self.nreads, self.outpath))
            #print args
#        print args
        result = [p.apply_async(window_core, arg) for arg in args]
        
#        result = p.map(window_core, args)
        #for arg in args:
        #    window_core(arg)
        #p = Process(target=self.window_core, args=(self.chrs_queue, self.nreads, self.outpath))
        #p.start()
        p.close()
        p.join()
        
def window_core(bamname, window_size, chr_tuple, nreads, outpath):
    bamfile = pysam.Samfile(bamname, 'rb')
    chr = chr_tuple[0]
    chr_length = chr_tuple[1]
    print chr
    wigfile = open(outpath + "/" + chr, 'w')
    out = "fixedStep chrom={0} start=1 step={1} span={1}\n".format(chr, window_size)
    wigfile.write(out)
    pos_vect = np.zeros((chr_length / window_size,))
    for read in bamfile.fetch(chr):
        read_mid = read.pos + read.isize / 2
        read_mid_index = read_mid / window_size
        pos_vect[read_mid_index] = pos_vect[read_mid_index] + 1
    for val in pos_vect:
        wigfile.write(str(10e6 * float(val) / nreads) + "\n")
    #out = "variableStep chrom={0}\n".format(chr)
    #wigfile.write(out)
    #curr_pos = 0
    #next_pos = curr_pos + window_size
    ##it = bamfile.pileup(chr)
    ##for column in it:
    ##    if column.pos % window_size == 0:
    ##        wigfile.write("\t".join([str(column.pos), str(10e6 * float(column.n) / nreads)]) + "\n")
    #count = 0        
    #while next_pos < chr_length:
    #    #print next_pos, chr_length
    #    count = bamfile.count(chr, curr_pos, next_pos)
    #    norm_count = 10e6 * float(count) / nreads
    #    print curr_pos, next_pos, chr_length, count, norm_count
    #    wigfile.write(str(norm_count))
    #    curr_pos = next_pos
    #    next_pos = curr_pos + window_size
#    count = 0

#    within = False
#    for read in bamfile.fetch(chr):
#        read_mid = read.pos + read.isize / 2
##        print read_mid, curr_pos, next_pos
#        while next_pos < chr_length:
#            if read_mid >= curr_pos and read_mid < next_pos:
#                print "hit"
#                count = count + 1
#                break
#            elif read_mid >= next_pos:
#                norm_count = count / nreads
#                wigfile.write(str(norm_count) + "\n")
#                count = 0
#            curr_pos = next_pos
#            next_pos = curr_pos + atoi(window_size)
#        print read_mid, curr_pos, next_pos, chr_length
#        sys.exit()
        """    
        print read_mid
        print curr_pos, next_pos
        if not within:
            curr_pos = next_pos + 1
            next_pos = curr_pos + atoi(window_size) - 1
        if read_mid > curr_pos and read_mid < next_pos:
            within = True
            count = count + 1
        elif read_mid > next_pos:
           
            curr_pos = next_pos + 1
            next_pos = curr_pos + atoi(window_size) - 1
            if next_pos >= chr_length: break
        """
        
    wigfile.close()

def main(argv):
    ## Take BAM of reads
    ## Count number of reads within specified window size
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', dest='bam')
    parser.add_argument('-w', dest='window')
    args = parser.parse_args()
    genome = open("/seq/lib/mouse.mm9.genome", 'r')
                 
    #bamfile = pysam.Samfile(args.bam, 'rb')
    #flagstat = open(args.bam + "_stat").readline()
    #nreads = flagstat.split()[0]
    outpath = "/".join([wigpath, os.path.basename(args.bam).split(".bam")[0], args.window])
    if not os.path.exists(outpath): os.makedirs(outpath)
    wi = windower(args.bam, outpath, args.window)
    wi.window()
    
if __name__ == "__main__":
    main(sys.argv)