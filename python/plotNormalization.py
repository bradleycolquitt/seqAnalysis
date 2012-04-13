#!/usr/bin/env python

import sys
import os
import argparse
import pdb
import tables as tb
import matplotlib.pyplot as mp
import itertools as it


        

class normCompare:
    def __init__(self, h5_track_name, cpg_track_name):
        self.h5_track = tb.openFile(h5_track_name)
        self.cpg_track = tb.openFile(cpg_track_name)
        #self.track_cals = 
        self.vals = {}
        self.n_bins = {}
        
    def processSamples(self, name):
        if name == "all":
            [self.collectValues(sample) for sample in self.h5_track.iterNodes("/")]
        else:
            self.collectValues(self.h5_track.getNode("/", name))
        
    def collectValues(self, sample):
        vals = self.vals
        n_bins = self.n_bins
        #pdb.set_trace()
        for sample_chr in sample._f_iterNodes():
            chr_name = sample_chr._v_name
            print chr_name
            cpg_chr = self.cpg_track.root.mm9._f_getChild(chr_name)
            for w, d in it.izip(cpg_chr, sample_chr):
                if w in vals:
                    vals[w] = vals[w] + d
                else:
                    vals[w] = d
                if w in n_bins:
                    n_bins[w] = n_bins[w] + 1
                else:
                    n_bins[w] = 1
    
    def plotCal(self):
        bins = self.vals.keys()
        cal_mean = []
        for bin in bins:
            cal_mean.append(self.vals[bin] / self.n_bins[bin])
        #cal_mean = [v / c for v,c in it.izip(self.vals, self.n_bins)]
        mp.scatter(bins, cal_mean)
        mp.show()
        
def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_data')
    parser.add_argument('--cpg_data')
    args = parser.parse_args()
    
    normObj = normCompare(args.sample_data, args.cpg_data)
    normObj.processSamples()
    return normObj

if __name__ == '__main__':
    main(sys.argv)