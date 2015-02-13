#!/usr/bin/env python

"""
Converts a BAM file to a WIG or Hdf5 track file by
    counting the number of fragments that lie within nonoverlapping genomic windows.
    A fragment is defined as a pair of reads for paired-end samples and
    an extended single-end read for single-end samples

Single-end samples: 
   Reads are extended by a given amount (param 'extend') and
   the midpoint of that extension is used to place a fragment in a genomic bin

Paired-end samples:
   The midpoint of a given pair of reads is used to determine fragment placement

Window counts are normalized by
    1. The number of reads (in millions) in the sample
    2. The window size to yield read number per kilobase
      - this allows comparison of values generated from different windows

WIGs and track files are outputted to separate directories
"""

import sys
import os
import re
import sam
import argparse
import pysam
import pdb
import file_util
import track_util
import signal_utils
import numpy as np
import numpy.ma as ma
import tables as tb
from subprocess import Popen
from scipy.stats.mstats import tmean
from multiprocessing import Queue, Process, Pool
from string import *

#bam_dir = "/media/storage2/data/rna/tophat"
#bam_dir = "/media/storage2/data/bam"
wig_dir = "/media/storage2/data/wig"
h5_dir = "/media/data2/hdf5"

class windower:
    def __init__(self, bamname, outname, track_name, 
                 window_size, extend, pe, pseudo,
                 ends, smooth, norm_by_mean, no_norm, output_type):
        
        self.bamname = bamname
        self.bamfile = pysam.AlignmentFile(bamname, 'rb')

        self.outname = outname
        self.output_type = output_type
        self.outfile = ""

        self.track_name = track_name
        if output_type == "wig":
            self.outfile = open(outname, 'w')
        elif output_type == "h5":
            self.outfile = tb.openFile(outname, "a")
            test = track_util.checkIfNodeExists(self.outfile, self.track_name, 
                                                create=True, accept_existence=False)
            if test: sys.exit()
        elif output_type == "none":
            pass

        self.window_size = atoi(window_size)
        self.pe = pe
        self.extend = int(extend)

        self.nreads = self.bamfile.mapped

        self.pseudo = pseudo
        self.ends = ends
        self.smooth = smooth
        self.norm_by_mean = norm_by_mean
        self.no_norm = no_norm

        #pdb.set_trace()
        self.window_correct = {}
        self.norm_path = ".".join([bamname, "chr_means_0"])
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
        
    # Interface for windowing functions
    def window(self):
        results = [window_core(self, chrom) for chrom in self.chrs_queue]
        #for result in results:
        #    self.window_correct[result[0]] = result[1]
     
    # Write window_correct dictionary if not already present
    def write_norm_vals(self):
        if not os.path.exists(self.norm_path):
            file = open(self.norm_path, 'w')
            for item in self.window_correct.iteritems():
                out = "\t".join([item[0], str(item[1])]) + "\n"
                file.write(out)

def window_core(obj, chrom):
    bamfile = pysam.AlignmentFile(obj.bamname, 'rb')
    chr_length = chrom[1]
    chrom = chrom[0]
    print chrom
        
    pos_vect = np.zeros((chr_length / obj.window_size,))
    
    read_mid = 0
    
    # Loop through BAM reads, adding midpoints to position vector
    for read in bamfile.fetch(chrom):
        read_mid = sam.read_mid_compute(obj, read)
        if read_mid < 0: continue
        read_mid_index = read_mid / obj.window_size
        if read_mid_index <= len(pos_vect) - 1:
            pos_vect[read_mid_index] = pos_vect[read_mid_index] + 1
    
    # Pseudocount 0 windows
    if obj.pseudo:
        pos_vect[pos_vect==0] = 1
    
    # Normalize by Reads per million and by reads per kilobase
    #window_correct = 1e6 * 1e3 / (float(obj.window_size) * obj.nreads)

    # Normalize by Reads per million
    window_correct = 1e6 * obj.nreads
    pos_vect = window_correct * pos_vect

    if obj.smooth > 0:
        pos_vect = signal_utils.smooth(pos_vect, obj.smooth, window="flat")

    if not obj.output_type == "none":
        write_values(obj, pos_vect, chrom, chr_length)

def write_values(obj, pos_vect, chrom, chr_length):
    if obj.output_type == "wig":
        out = "fixedStep chrom={0} start={1} step={2} span={2}\n".format(chrom, 
                                                                         "1", 
                                                                         obj.window_size)
        obj.outfile.write(out)
        for val in pos_vect:
            obj.outfile.write(str(val) + "\n")
    elif obj.output_type == "h5":
        obj.outfile.createArray("/" + obj.track_name, chrom, pos_vect)
        track_util.setTrackAttributes(obj.outfile, 
                                      obj.track_name + "/" + chrom, 
                                      0, 
                                      chr_length / obj.window_size, 
                                      chrom, 
                                      obj.window_size)

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', dest='bam', required=False)
    parser.add_argument('-w', dest='window')
    parser.add_argument('-e', dest='extend', type=int, required=False, default=0)
    parser.add_argument('-t', dest='output_type', default="wig", choices=['wig', 'h5', 'none'])
    parser.add_argument('--h5_filename')
    parser.add_argument('--paired_end', action='store_true', default=False, help="Is sample PE? If set with 'extend', will extend by given amount from read1 only.")
    parser.add_argument('--pseudocount', dest='pseudo', action='store_true', default=False)
    parser.add_argument('--ends', action='store_true', default=False, help="Record position of the end of each read")
    parser.add_argument('--smooth', type=int, help="Size of smoothing window")
    parser.add_argument('--norm_by_mean', action='store_true', default=False, help="Normalize values by mean across chromosome (excluding zero windows), not total read number")
    parser.add_argument('--no_norm', action='store_true', default=False, help="Don't normalize. Just count number of reads.")
    args = parser.parse_args()

    #bam_path = "/".join([bam_dir, args.bam])
    bam_prefix = os.path.basename(args.bam).split(".bam")[0]

    out_dir = wig_dir
    if args.output_type == "h5":
        out_dir = h5_dir
    if args.norm_by_mean:
       out_path = "/".join([out_dir, "mean"])
    elif args.no_norm:
        out_path = "/".join([out_dir, "no_norm"])
    else:
        out_path = "/".join([out_dir, "rpm"])

    track_name = ""
    outfile = ""
    if args.output_type == "wig":
        out_path = "/".join([out_path, bam_prefix])
        outfile = "/".join([out_path, "_".join([bam_prefix, args.window])])

        if args.extend > 0:
            outfile = "_".join([outfile, "extend" + str(args.extend)])
        if args.pseudo:
            outfile = "_".join([outfile, "pseudo"])
        if args.ends:
            outfile = "_".join([outfile, "ends"])
        if args.smooth:
            outfile = "_".join([outfile, "smooth" + str(args.smooth)])
        outfile = outfile + ".wig"
        print outfile
        if not args.output_type == "none":
            if not os.path.exists(out_path): os.makedirs(out_path)
            if os.path.exists(outfile):
                dec = raw_input("WIG exists. Overwrite [y/n]? ")
                if dec == "n": return

    elif args.output_type == "h5":
        outfile = "/".join([out_path, "_".join([args.h5_filename, args.window])]) + ".h5"
        print outfile
        track_name = bam_prefix
        if args.extend > 0:
            track_name = "_".join([track_name, "extend" + str(args.extend)])
        if args.pseudo:
            track_name = "_".join([track_name, "pseudo"])
        if args.ends:
            track_name = "_".join([track_name, "ends"])
        if args.smooth:
            track_name = "_".join([track_name, "smooth" + str(args.smooth)])

    wi = windower(args.bam, outfile, track_name,
                  args.window, args.extend, args.paired_end,
                  args.pseudo, args.ends, args.smooth,
                  args.norm_by_mean, args.no_norm, args.output_type)
    wi.window()

    if args.output_type == "h5":
        wi.outfile.flush()

if __name__ == "__main__":
    main(sys.argv)
