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
import signal_utils
import numpy as np
from subprocess import Popen
from multiprocessing import Queue, Process, Pool
from string import *

#bam_dir = "/home/user/data/rna/tophat"
bam_dir = "/media/storage2/data/bam"
#wig_dir = "/media/storage2/data/wig/unnorm"

#bam_dir = "/home/user/storage/data/guanz/bam"
wig_dir = "/media/storage2/data/wig"
tdf_dir = "/media/storage2/data/tdf"

class windower:
    def __init__(self, bamname, wigname, window_size, extend, pe, bed, pseudo,
                 full, ends, smooth, norm_by_mean, no_norm, no_output):
        
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
        self.smooth = smooth
        self.norm_by_mean = norm_by_mean
        self.no_output = no_output
        
        #pdb.set_trace()
        self.window_correct = {}
        self.norm_path = ".".join([bamname, "chr_means"])
        if os.path.exists(self.norm_path):
            for line in open(self.norm_path, 'r'):
                sline = line.split()
                self.window_correct[sline[0]] = float(sline[1])
        
        if pe:
            # If strand separated BAM
            if re.search("plus", bamname) or re.search("minus", bamname):
                self.nreads = self.bamfile.mapped
            # Normalize to number of fragments
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
        results = 0
        if self.bed_name:
                window_bed(self)      
        else:
            if not self.full:
                results = [window_core(self, chrom) for chrom in self.chrs_queue]
            else:
                results = [window_full(self, chrom) for chrom in self.chrs_queue]
            
            #for chrom in self.chrs_queue:
            #    if not self.full:
            #        window_core(self, chrom)
            #    else:
            #        window_full(self, chrom)
        #pdb.set_trace()
        for result in results:
            self.window_correct[result[0]] = result[1]
        """    
        args = []
        for chr in self.chrs_queue:
            args.append((self.bamname, self.wigfile, self.window_size, chr,
                         self.pe, self.extend, self.nreads, self.bed_file,
                         self.pseudo, self.ends, self.smooth, self.norm_by_mean))
        for arg in args:
            
            if not self.full:
                #getattr(self, window_core(args))
                #window_core(arg[0], arg[1], arg[2], arg[3], arg[4], arg[5],
                #            arg[6], arg[7], arg[8], arg[9], arg[10])
                window_core(arg)
            else:
                #window_full(arg[0], arg[1], arg[2], arg[3], arg[4], arg[5],
                #            arg[6], arg[7], arg[)
                window_full(self)
        self.wigfile.close()
        """
    # Write window_correct dictionary if not already present
    def write_norm_vals(self):
        if not os.path.exists(self.norm_path):
            file = open(self.norm_path, 'w')
            for item in self.window_correct.iteritems():
                out = "\t".join([item[0], str(item[1])]) + "\n"
                file.write(out)
    
    # Convert WIG to TDF            
    def tdf(self):
        print self.tdffile
        cmd_args = ['igvtools', 'tile', self.wigname, self.tdffile, 'mm9']
        p = Popen(cmd_args)
        p.wait()
        
#def window_core(bamname, wigfile, window_size, chr_tuple, pe, extend, nreads,
#                bed, pseudo, ends, smooth):
def window_core(obj, chrom):    
    
    """
    bamname = arg[0]
    wigfile = arg[1]
    window_size = arg[2]
    chr_tuple = arg[3]
    pe = arg[4]
    extend = arg[5]
    nreads = arg[6]
    bed = arg[7]
    pseudo = arg[8]
    ends = arg[9]
    smooth = arg[10]
    norm_by_mean = arg[11]
    """
    
    bamfile = pysam.Samfile(obj.bamname, 'rb')
    #pdb.set_trace()
    chr_length = chrom[1]
    chrom = chrom[0]
    print chrom
    
    # Normalize by Reads per million and by reads per kilobase
    window_correct = 1e6 * 1e3 / (float(obj.window_size) * obj.nreads)
        
    pos_vect = np.zeros((chr_length / obj.window_size,))
    
    read_mid = 0
    
    for read in bamfile.fetch(chrom):
        
        # if just using read ends
        if obj.ends:
           # pdb.set_trace()
            if read.is_reverse:
                read_mid = read.aend + 1
            else:
                read_mid = read.pos + 2
        elif obj.pe:
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
        read_mid_index = read_mid / obj.window_size
        if read_mid_index <= len(pos_vect) - 1:
            pos_vect[read_mid_index] = pos_vect[read_mid_index] + 1
    
    #pdb.set_trace()
    first_non_zero = np.nonzero(pos_vect)[0][0]
    first_non_zero_scaled = first_non_zero * obj.window_size
    
    out = "fixedStep chrom={0} start={1} step={2} span={2}\n".format(chrom, first_non_zero_scaled, obj.window_size)
    obj.wigfile.write(out)
    
    pos_vect = window_correct * pos_vect
    #pdb.set_trace()
    if obj.smooth > 0:
        pos_vect = signal_utils.smooth(pos_vect, obj.smooth, window="flat")
    
    if not obj.no_norm:    
        for val in pos_vect[first_non_zero:]:
            if obj.pseudo and val == 0: val = 1
            obj.wigfile.write(str(val) + "\n")
        
def window_full(obj, chrom):
    #pdb.set_trace()
   
    bamfile = pysam.Samfile(obj.bamname, 'rb')
    chr_length = chrom[1]
    chrom = chrom[0]
    print chrom
    
    window_correct = 0
    if not obj.norm_by_mean:
        window_correct = (1e6) / (float(obj.window_size) * obj.nreads)
    out = "fixedStep chrom={0} start=1 step={1} span={1}\n".format(chrom, obj.window_size)
    obj.wigfile.write(out)
    
    ## Initialize output vector and BAM iterator
    pos_vect_len = chr_length / obj.window_size
    pos_vect = np.zeros((pos_vect_len,))
#    pos_vect = np.zeros((obj.chr_lengths[chrom] / obj.window_size,))
    #it = bamfile.pileup(reference=chr)
    start = 0
    end = 0
    
    ## Loop through each position in iterator
    for read in bamfile.fetch(chrom):
#        pdb.set_trace()
        if obj.pe:
            if not read.is_reverse: 
                start = read.pos
                end = read.pnext + read.qlen
                #print "{0},{1},{2}".format(start, end, (end-start))
                if (end - start) > 2000:
                    # pdb.set_trace()
                    continue
            else:
                continue
        if obj.extend:
            if not read.is_reverse:
                start = read.pos
                end = start + extend
            else:
                end = read.aend
                start = end - extend
        
        start /= obj.window_size
        end /= obj.window_size
        
        if end >= pos_vect_len: end = pos_vect_len - 1        
        try:
            for i in range(start, end): pos_vect[i] += 1
        except IndexError:
            pdb.set_trace()
    
    
    if obj.norm_by_mean:
        mean = np.mean(pos_vect)
        if mean > 0:
            window_correct = 1 / float(mean)
    if not obj.no_output:
        print "Start write"
        for val in pos_vect:
            if obj.pseudo and val == 0: val = 1
            obj.wigfile.write(str(window_correct * float(val)) + "\n")
    
    return((chrom, window_correct))
    
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
        
        bed_start = int(line[1]) - 1
        bed_end = int(line[2]) - 1
        
        # Setup position vector
        pos_vect = np.zeros(((bed_end - bed_start) / obj.window_size,))

        # Loop through intersecting reads
        if not obj.full:
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
        else:
            # Currently written for only window_size of 1
#            pdb.set_trace()
            for column in bamfile.pileup(reference=chrom, start=bed_start, end=bed_end):
                if (column.pos >= bed_start and column.pos < bed_end):
 #                   pdb.set_trace()
                    try:
                        pos_vect[(column.pos - bed_start)] = column.n
                    except:
                        pdb.set_trace()
                        
        if obj.norm_by_mean:
            if not chrom in obj.window_correct:
                total_n = 0
                count_n = 0
                print "Computing average values..."
                for column in bamfile.pileup(reference=chrom):
                    total_n += column.n
                    count_n += 1
                
                obj.window_correct[chrom] = 1 / (total_n / float(count_n))
                print "Average coverage = {0}".format(obj.window_correct[chrom])
            window_correct = obj.window_correct[chrom]
            
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
    parser.add_argument('--smooth', type=int, help="Size of smoothing window")
    parser.add_argument('--norm_by_mean', action='store_true', default=False, help="Normalize values by mean across chromosome (excluding zero windows), not total read number")
    parser.add_argument('--no_norm', action='store_true', default=False, help="Don't normalize. Just count number of reads.")
    parser.add_argument('--no_output', action='store_true', default=False, help="Don't output WIG. Useful if just computing norm values")
    args = parser.parse_args()
      
    bam_path = "/".join([bam_dir, args.bam])
    bam_prefix = os.path.basename(args.bam).split(".bam")[0]
    
    if args.norm_by_mean:
       wig_path = "/".join([wig_dir, "mean"])
    elif args.no_norm:
        wig_path = "/".join([wig_dir, "no_norm"])
    else:
        wig_path = "/".join([wig_dir, "rpkm"])
        
    wig_path = "/".join([wig_path, bam_prefix])
    wig_file = "/".join([wig_path, "_".join([bam_prefix, args.window])])

    if args.bed:
        wig_file = "_".join([wig_file, os.path.basename(args.bed)])
    if args.pseudo:
        wig_file = "_".join([wig_file, "pseudo"])
    if args.full:
        wig_file = "_".join([wig_file, "full"])
    if args.ends:
        wig_file = "_".join([wig_file, "ends"])
    if args.smooth:
        wig_file = "_".join([wig_file, "smooth" + str(args.smooth)])
    print wig_file
    wig_file = wig_file + ".wig"

    if not os.path.exists(wig_path): os.makedirs(wig_path)
    if os.path.exists(wig_file):
        dec = raw_input("WIG exists. Overwrite [y/n]? ")
        if dec == "n": return

    wi = windower(bam_path, wig_file, args.window, args.extend, args.paired_end,
                  args.bed, args.pseudo, args.full, args.ends, args.smooth,
                  args.norm_by_mean, args.no_norm, args.no_output)
    wi.window()
    wi.write_norm_vals()
    wi.wigfile.close()
            
if __name__ == "__main__":
    main(sys.argv)
