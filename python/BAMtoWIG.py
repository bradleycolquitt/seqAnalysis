#!/usr/bin/env python

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
bam_dir = "/media/storage2/data/bam"
#bam_dir = "/media/storage2/data/rna/tophat"
#wig_dir = "/media/storage2/data/wig/unnorm"
wig_dir = "/media/storage2/data/wig/rpkm"
tdf_dir = "/media/storage2/data/tdf"

class windower:
    def __init__(self, bamname, wigname, window_size, extend, pe, pseudo, full):
        #print pe
        self.bamname = bamname
        self.bamfile = pysam.Samfile(bamname, 'rb')
        self.wigname = wigname
        self.wigfile = open(wigname, 'w')
        self.tdffile = "/".join([tdf_dir, "".join([os.path.basename(wigname).split(".wig")[0], ".tdf"])])
        self.window_size = atoi(window_size)
        self.pe = pe
        self.extend = int(extend)
        #self.nreads = 0
        self.nreads = self.bamfile.mapped
        self.pseudo = pseudo
        self.full = full
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
        p = Pool(processes=6)
        args = []
        for chr in self.chrs_queue:
            args.append((self.bamname, self.wigfile, self.window_size, chr, self.pe, self.extend, self.nreads, self.pseudo))
        for arg in args:
            if not self.full:
                window_core(arg[0], arg[1], arg[2], arg[3], arg[4], arg[5], arg[6], arg[7])
            else:
                window_full(arg[0], arg[1], arg[2], arg[3], arg[4], arg[5], arg[6], arg[7])
#                result = [p.apply_async(window_full, arg) for arg in args]        
        #print args
        
#        p.close()
#        p.join()
        #pool = Pool(processes=5)
        #bam = pysam.Samfile(self.sample, 'rb')
        #for chr_tbp in chrs_tbp:
        #    pool.apply_async(tab_bam, (self.bamname, self.wigfile, self.window_size, chr, self.pe, self.extend, self.nreads, self.pseudo))))
            #tab_bam(self, chr_tbp, tmp_path)
        #pool.close()
        #pool.join()
        self.wigfile.close()
        
    def tdf(self):
        print self.tdffile
        cmd_args = ['igvtools', 'tile', self.wigname, self.tdffile, 'mm9']
        p = Popen(cmd_args)
        p.wait()
        
def window_core(bamname, wigfile, window_size, chr_tuple, pe, extend, nreads, pseudo):
    bamfile = pysam.Samfile(bamname, 'rb')
    chr = chr_tuple[0]
    chr_length = chr_tuple[1]
    print chr
    window_correct = 1000 / float(window_size)
    #wigfile = open(outpath + "/" + chr, 'w')
    out = "fixedStep chrom={0} start=1 step={1} span={1}\n".format(chr, window_size)
    wigfile.write(out)
    
    pos_vect = np.zeros((chr_length / window_size,))
    #print pe
    #print extend
    for read in bamfile.fetch(chr):
        #print read.pos
        read_mid = 0
        if pe:
            if read.is_read1:
                read_mid = read.pos + read.isize / 2
            else: continue
        else:
            read_mid = read.pos + extend / 2
            #print read_mid
        read_mid_index = read_mid / window_size
        if read_mid_index <= len(pos_vect) - 1:
            pos_vect[read_mid_index] = pos_vect[read_mid_index] + 1
    for val in pos_vect:
        if pseudo and val == 0: val = 1
        #wigfile.write(str((window_correct * 1e6 * float(val)) / nreads) + "\n")
        wigfile.write(str(window_correct * float(val)) + "\n")
       
    #wigfile.close()
    
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
    #for proxy in it:
    #    read_index = proxy.pos / window_size
    #    if read_index <= len(pos_vect) - 1:
    #        pos_vect[read_index] = pos_vect[read_index] + proxy.n
    print "Start loop"
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
        
       
def main(argv):
    ## Take BAM
    ## Count number of reads within specified window size
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', dest='bam', required=False)
    parser.add_argument('-w', dest='window')
    
    parser.add_argument('-e', dest='extend', type=int, required=False, default=0)
    parser.add_argument('-d', dest='date', required=False)
    parser.add_argument('--paired_end', action='store_true', default=False)
    parser.add_argument('--pseudocount', dest='pseudo', action='store_true', default=False)
    parser.add_argument('--full', action='store_true', default=False, help="Record extent of each read")
    args = parser.parse_args()
    
        
        
    if args.date:
        bams = [f for f in os.listdir("/".join([bam_dir, args.date])) if re.search(".bam$", f)]
        for bam in bams:
            bam_path = "/".join([bam_dir, args.date, bam])
            print bam
            bam_prefix = os.path.basename(bam).split(".bam")[0]
            wig_path = "/".join([wig_dir, bam_prefix])
            wig_file = ""
            if not args.pseudo:
                wig_file = "/".join([wig_path, "_".join([bam_prefix, args.window]) + ".wig"])
            else:
                wig_file = "/".join([wig_path, "_".join([bam_prefix, args.window]) + "pseudo.wig"])
            if not os.path.exists(wig_path): os.makedirs(wig_path)
            if os.path.exists(wig_file):
                dec = raw_input("WIG exists. Overwrite [y/n]? ")
                if dec == "n": continue
            
            wi = windower(bam_path, wig_file, args.window, args.extend, args.paired_end, args.pseudo, args.full)
            wi.window()
            wi.wigfile.close()
            #wi.tdf()
    else:
        bam_path = "/".join([bam_dir, args.bam])
        bam_prefix = os.path.basename(args.bam).split(".bam")[0]
        wig_path = "/".join([wig_dir, bam_prefix])
        wig_file = "/".join([wig_path, "_".join([bam_prefix, args.window])])
        
        if args.pseudo:
            wig_file = "_".join([wig_file, "pseudo"])
        if args.full:
            wig_file = "_".join([wig_file, "full"])
            
        wig_file = wig_file + ".wig"
        
        if not os.path.exists(wig_path): os.makedirs(wig_path)
        if os.path.exists(wig_file):
            dec = raw_input("WIG exists. Overwrite [y/n]? ")
            if dec == "n": return
            
        wi = windower(bam_path, wig_file, args.window, args.extend, args.paired_end, args.pseudo, args.full)
        wi.window()
        wi.wigfile.close()
        #wi.tdf()
        
            
    
if __name__ == "__main__":
    main(sys.argv)