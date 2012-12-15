#!/usr/bin/env python

"""
Converts a BAM file to a WIG file by
    counting the number of fragments that lie within nonoverlapping genomic windows.
    A fragment is defined as a pair of reads for paired-end samples and
    an extended single-end read for single-end samles

Single-end samples: 
   Reads are extended by a given amount (param 'extend') and
   the midpoint of that extension is used to place a fragment in a genomic bin

Paired-end samples:
   The midpoint of a given pair of reads is used to determine fragment placement

Window counts are normalized by
    1. The number of reads (in millions) in the sample
    2. The window size to yield read number per kilobase
      - this allows comparison of values generated from different windows

WIGs are outputted to a common directory
"""

import sys
import os
import re
import argparse
import pysam
import pdb
import file_util
import numpy as np
from subprocess import Popen
from multiprocessing import Queue, Process, Pool
from string import *

#bam_dir = "/home/user/data/rna/tophat"
#bam_dir = "/media/storage2/data/bam"
bam_dir = "/media/storage2/data/rna/tophat"
#wig_dir = "/media/storage2/data/wig/unnorm"
wig_dir = "/media/storage2/data/wig/rpkm"
tdf_dir = "/media/storage2/data/tdf"

class windower:
    def __init__(self, bamname, wigname, window_size, extend, pe, bed, pseudo, full, ends):
        #print pe
        self.bamname = bamname
        self.bamfile = pysam.Samfile(bamname, 'rb')
        self.wigname = wigname
        self.wigfile = open(wigname, 'w')
        self.tdffile = "/".join([tdf_dir, "".join([os.path.basename(wigname).split(".wig")[0], ".tdf"])])
        self.window_size = atoi(window_size)
        self.pe = pe
        self.extend = int(extend)

        self.nreads = self.bamfile.mapped
        self.bed_name = bed
        self.bed_file = ""
        if bed:
            self.bed_file = open(bed)
        self.pseudo = pseudo
        self.full = full
        self.ends = ends
        
        if pe:
            if re.search("plus", bamname) or re.search("minus", bamname):
                self.nreads = self.bamfile.mapped
            else:    
                self.nreads = self.bamfile.mapped / 2
        else:
            self.nreads = self.bamfile.mapped
        print self.nreads
        self.chr_lengths = self.bamfile.lengths
        self.chrs_queue = []
        for index in range(self.bamfile.nreferences):
            self.chrs_queue.append((self.bamfile.references[index], self.bamfile.lengths[index]))
        
    def window(self):
        if self.bed_name:
                window_bed(self)
                return
        #p = Pool(processes=6)
        args = []
        for chr in self.chrs_queue:
            args.append((self.bamname, self.wigfile, self.window_size, chr, self.pe,
                         self.extend, self.nreads, self.bed_file, self.pseudo, self.ends))
        for arg in args:
            
            if not self.full:
                window_core(arg[0], arg[1], arg[2], arg[3], arg[4], arg[5], arg[6], arg[7], arg[8], arg[9])
            else:
                window_full(arg[0], arg[1], arg[2], arg[3], arg[4], arg[5], arg[6], arg[7])
        
        self.wigfile.close()
        
    def tdf(self):
        print self.tdffile
        cmd_args = ['igvtools', 'tile', self.wigname, self.tdffile, 'mm9']
        p = Popen(cmd_args)
        p.wait()
        
def window_core(bamname, wigfile, window_size, chr_tuple, pe, extend, nreads, bed, pseudo, ends):
    bamfile = pysam.Samfile(bamname, 'rb')
    chr = chr_tuple[0]
    chr_length = chr_tuple[1]
    print chr
    #window_correct = 1000 / float(window_size)
    
    # Normalize by Reads per million and by reads per kilobase
    window_correct = 1e6 * 1e3 / (float(window_size) * nreads)
    
    
    pos_vect = np.zeros((chr_length / window_size,))
    
    read_mid = 0
    
    for read in bamfile.fetch(chr):
        #read_mid = 0
        
        # if just using read ends
        if ends:
            if read.is_reverse:
                read_mid = read.aend - 1
            else:
                read_mid = read.pos
        elif pe and not ends:
            if read.is_read1:
                if read.is_reverse:
                    read_mid = read.aend + read.isize / 2
                else:
                    read_mid = read.pos + read.isize / 2
            else: continue
        else:
            if read.is_reverse:
                read_mid = read.aend - extend / 2
            else:
                read_mid = read.pos + extend / 2
        read_mid_index = read_mid / window_size
        if read_mid_index <= len(pos_vect) - 1:
            pos_vect[read_mid_index] = pos_vect[read_mid_index] + 1
    
    #pdb.set_trace()
    first_non_zero = np.nonzero(pos_vect)[0][0]
    first_non_zero_scaled = first_non_zero * window_size
    
    out = "fixedStep chrom={0} start={1} step={2} span={2}\n".format(chr, first_non_zero_scaled, window_size)
    wigfile.write(out)
    
    for val in pos_vect[first_non_zero:]:
        if pseudo and val == 0: val = 1
        wigfile.write(str(window_correct * float(val)) + "\n")
        
def window_full(bamname, wigfile, window_size, chr_tuple, pe, extend, nreads, pseudo):
    #pdb.set_trace()
    bamfile = pysam.Samfile(bamname, 'rb')
    chr = chr_tuple[0]
    chr_length = chr_tuple[1]
    print chr
    #window_correct = 1000 / float(window_size)
    window_correct = (1e6) / (float(window_size) * nreads)
    out = "fixedStep chrom={0} start=1 step={1} span={1}\n".format(chr, window_size)
    wigfile.write(out)
    
    ## Initialize output vector and BAM iterator
    pos_vect = np.zeros((chr_length / window_size,))
    #it = bamfile.pileup(reference=chr)
    start = 0
    end = 0
    
    ## Loop through each position in iterator
    for read in bamfile.fetch(chr):
#        pdb.set_trace()
        if pe:
            if not read.is_reverse:
                start = read.pos
                end = read.pnext + read.qlen
        if extend:
            if not read.is_reverse:
                start = read.pos
                end = start + extend
            else:
                end = read.aend
                start = end - extend
        start /= window_size
        end /= window_size
        if end >= len(pos_vect): end = len(pos_vect) - 1
        
        try:
            for i in range(start, end): pos_vect[i] += 1
        except IndexError:
            pdb.set_trace()
    print "Start write"    
    for val in pos_vect:
        if pseudo and val == 0: val = 1
        wigfile.write(str(window_correct * float(val)) + "\n")
        
def window_bed(obj):
    
    bamfile = pysam.Samfile(obj.bamname, 'rb')
    #chr = chr_info[0]
    #chr_length = chr_info[1]
    #print chr
        
    # Normalize by Reads per million and by reads per kilobase
    window_correct = 1e6 * 1e3 / (float(obj.window_size) * obj.nreads)
    
    line = ""
    chrom = ""
    curr_chrom = ""
    bed_start = 0
    bed_end = 0
    pos_vect = ""
    read_mid = 0
    read_mid_index = 0
    
    # Loop through bed regions
    for line in obj.bed_file:
        line = line.split()
        chrom = line[0]
        if chrom != curr_chrom:
            print chrom
            curr_chrom = chrom
        
        bed_start = int(line[1])
        bed_end = int(line[2])
        
        # Setup position vector
        pos_vect = np.zeros(((bed_end - bed_start) / obj.window_size,))

        # Loop through intersecting reads
        for read in bamfile.fetch(reference=chrom, start=bed_start, end=bed_end):
            # Extract ends
            if obj.ends:
                if read.is_reverse:
                    read_mid = read.aend - 1 - bed_start
                else:
                    read_mid = read.pos - bed_start
            
            # Record ends within relative vector
            read_mid_index = read_mid / obj.window_size
            if read_mid_index >= 0 and read_mid_index <= len(pos_vect) - 1:
                pos_vect[read_mid_index] = pos_vect[read_mid_index] + 1
    
        # Write WIG, header for each bed region
        out = "fixedStep chrom={0} start={1} step={2} span={2}\n".format(chrom, str(bed_start), obj.window_size)
        obj.wigfile.write(out)
    
        for val in pos_vect:
            obj.wigfile.write(str(window_correct * float(val)) + "\n")

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', dest='bam', required=False)
    parser.add_argument('-w', dest='window')
    parser.add_argument('-e', dest='extend', type=int, required=False, default=0)
    parser.add_argument('--bed')
    parser.add_argument('--paired_end', action='store_true', default=False)
    parser.add_argument('--pseudocount', dest='pseudo', action='store_true', default=False)
    parser.add_argument('--full', action='store_true', default=False, help="Record extent of each read")
    parser.add_argument('--ends', action='store_true', default=False, help="Record position of the end of each read")
    
    args = parser.parse_args()
      
    bam_path = "/".join([bam_dir, args.bam])
    bam_prefix = os.path.basename(args.bam).split(".bam")[0]
    wig_path = "/".join([wig_dir, bam_prefix])
    wig_file = "/".join([wig_path, "_".join([bam_prefix, args.window])])

    if args.bed:
        wig_file = "_".join([wig_file, os.path.basename(args.bed)])
    if args.pseudo:
        wig_file = "_".join([wig_file, "pseudo"])
    if args.full:
        wig_file = "_".join([wig_file, "full"])
    if args.ends:
        wig_file = "_".join([wig_file, "ends"])
    print wig_file
    wig_file = wig_file + ".wig"

    if not os.path.exists(wig_path): os.makedirs(wig_path)
    if os.path.exists(wig_file):
        dec = raw_input("WIG exists. Overwrite [y/n]? ")
        if dec == "n": return

    wi = windower(bam_path, wig_file, args.window, args.extend, args.paired_end, args.bed, args.pseudo, args.full, args.ends)
    wi.window()
    wi.wigfile.close()
            
if __name__ == "__main__":
    main(sys.argv)
