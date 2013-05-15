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
import itertools
import pysam
import signal_utils
import math
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

        # Remove sex chromosomes to account for genotype differences
        chr_tbp = [a._v_name for a in sample1._f_listNodes()]
        chr_tbp = filter(lambda x: not re.search("X|Y|M", x), chr_tbp)
        
        for chrom in chr_tbp:
            chrom_data = sample1._f_getChild(chrom)[:]
            sample1_data.append(chrom_data)

            chrom_data = sample2._f_getChild(chrom)[:]
            sample2_data.append(chrom_data)
        
        # Concatenate chromosome data into one vector    
        vals1 = np.concatenate(sample1_data)
        vals2 = np.concatenate(sample2_data)
        
        # Reduce resolution by specified number of windows
        if self.res_reduce > 0:
            column = len(vals1) / self.res_reduce
            total = self.res_reduce * column
            vals1 = vals1[:total].reshape(column, self.res_reduce).mean(axis=1)
            vals2 = vals2[:total].reshape(column, self.res_reduce).mean(axis=1)
        
        # Group into groups of 100 (array of Nx100)
        row = len(vals1) / 100
        total = 100 * row
        vals1_mat = vals1[:total].reshape(len(vals1) / 100, 100)
        vals2_mat = vals2[:total].reshape(len(vals2) / 100, 100)
        corr = 0

        # Accumulate correlation values for each group
        for r1, r2 in itertools.izip(vals1_mat, vals2_mat):
            if self.corr_type == "pearsonr":
                val = stats.pearsonr(r1, r2)[0]
            elif self.corr_type == "spearmanr":
                val = stats.spearmanr(r1, r2)[0]
            if math.isnan(val): continue
            corr = corr + val

        corr = corr / (len(vals1) / 100)        
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
