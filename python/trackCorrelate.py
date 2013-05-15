#! /usr/bin/env python

####
# Correlate (cross for two samples auto for one) values across given feature
# Outputs
####
import sys
import os
import re
import shutil
import argparse
import tempfile
import pdb
import pysam
import signal_utils
import numpy as np
import scipy.signal as ss
import scipy.stats as stats
import tables as tb
from multiprocessing import Pool

feature_path = "/home/user/lib/features_general/"
out_path = "/media/storage2/analysis/crosscor"

class cross_track:
    
    def __init__(self, file_a, file_b, track_a, track_b, res_reduce, corr_type):
        self.file_a = file_a
        self.file_b = file_b
        self.track_a = track_a
        self.track_b = track_b
        self.res_reduce = res_reduce
        self.corr_type = corr_type    

    def run(self):

        h5_a = tb.openFile(self.file_a)
        h5_b = tb.openFile(self.file_b)
        sample1 = h5_a.getNode("/", self.track_a)
        sample_chrs = [chr._v_name for chr in sample1._f_iterNodes()]
        sample2 = h5_b.getNode("/", self.track_b)
        sample1_data = []
        sample2_data = []

        chr_tbp = [a._v_name for a in sample1._f_listNodes()]
        chr_tbp = filter(lambda x: not re.search("X|Y|M", x), chr_tbp)
        
        for chrom in chr_tbp:
            chrom_data = sample1._f_getChild(chrom)[:]
            chrom_data = signal_utils.smooth(chrom_data, self.res_reduce, window="flat")
            sample1_data.append(chrom_data)

            chrom_data = sample2._f_getChild(chrom)[:]
            chrom_data = signal_utils.smooth(chrom_data, self.res_reduce, window="flat")
            sample2_data.append(chrom_data)
        
        vals1 = np.concatenate(sample1_data)
        vals2 = np.concatenate(sample2_data)
        
        #pdb.set_trace()   
        cutoff = 1E3
        vals1[vals1>cutoff] = cutoff
        vals2[vals2>cutoff] = cutoff
        
        #vals1 = signal_utils.smooth(vals1, self.res_reduce, window="flat")
        #vals2 = signal_utils.smooth(vals2, self.res_reduce, window="flat")
        """
        if self.res_reduce > 0:
            column = len(vals1) / self.res_reduce
            total = self.res_reduce * column
            vals1 = vals1[:total].reshape(self.res_reduce, column).mean(axis=0)
            vals2 = vals2[:total].reshape(self.res_reduce, column).mean(axis=0)
        """
        not_zero = np.multiply(vals1, vals2) > 0
        
        if self.corr_type == "cross" or self.corr_type == "auto": 
            corr = ss.fftconvolve(vals1, vals2, 'same')
        elif self.corr_type == "spearmanr":
            corr = stats.spearmanr(vals1[not_zero], vals2[not_zero])
        elif self.corr_type == "pearsonr":
            print "here"
            corr = stats.pearsonr(vals1[not_zero], vals2[not_zero])

        print corr
    
def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--file_a", dest="file_a", required=False)
    parser.add_argument("--track_a", dest="track_a", required=False)
    parser.add_argument("--file_b", required=False)
    parser.add_argument("--track_b", dest="track_b", required=False)
    parser.add_argument("--reduce", type=int, default=0)
    parser.add_argument("--cor_type", choices=['cross', 'spearmanr', 'pearsonr'])
    
    args = parser.parse_args()
    
    file_b = args.file_a
    if args.file_b: file_b = args.file_b
    
    track_b = ""
    if args.track_b or args.file_b:
        # cross-correlate
        track_b = args.track_b
    else:
        # auto-correlate
        track_b = args.track_a
        
    obj = cross_track(args.file_a, file_b, args.track_a, track_b, args.reduce, args.cor_type)
    obj.run()

    

if __name__ == "__main__":
    main(sys.argv)
