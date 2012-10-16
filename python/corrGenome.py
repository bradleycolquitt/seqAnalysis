#! /usr/bin/env python

import sys
import os
import shutil
import argparse
import tempfile
import pdb
import math
import numpy as np
import scipy.stats.stats as ss
import tables as tb
from multiprocessing import Pool

feature_path = "/home/user/lib/features_general/"
out_path = "/media/storage2/analysis/genomecor"

## Correlate two tracks in sliding windows across genome
## Input h5, two tracks, window size, step size
## Output as wig

class corr_genome:
    def __init__(self, h5, track1, track2, data_type, method, window, step):
        self.h5 = h5
        self.track1 = track1
        self.track2 = track2
        self.func = ss.pearsonr
        if method == "spearman": self.func = ss.spearmanr
        elif method == "ratio": self.func = val_ratio
        elif method == "diff": self.func = val_diff
        elif method == "fraction": self.func = val_fraction
        self.window = int(window)
        self.step = int(step)
        self.tmp_path = tempfile.mkdtemp()
        self.out_path = "/".join([out_path, data_type, method])
        if not os.path.exists(self.out_path): os.makedirs(self.out_path)
        self.outfile = "/".join([self.out_path, "_".join([track1, track2, "W" + window, "S" + step]) + ".wig"])
        if os.path.exists(self.outfile):
            dec = raw_input("File exists. Overwrite? [y/n]")
            if dec == "n": sys.exit()
        
        # get chrs tbp
        h5_data = tb.openFile(h5)
        sample1 = h5_data.getNode("/", track1)
        sample2 = h5_data.getNode("/", track2)
        sample1_chrs = [chr._v_name for chr in sample1._f_iterNodes()]
        sample2_chrs = [chr._v_name for chr in sample2._f_iterNodes()]
        #pdb.set_trace()
        self.chrs_tbp = list(set(sample1_chrs) & set(sample2_chrs))
        
        #if not chr_tbp in sample_chrs: return
        
        
    def file_combine(self):        
        out = open(self.outfile, 'w')
        files_tmp = os.listdir(self.tmp_path)
        for file_tmp in files_tmp:
            a = open(self.tmp_path + "/" + file_tmp)
            for line in a:
                out.write(line)
        shutil.rmtree(self.tmp_path)
        out.close()

def compute_by_chrom(obj):
    
        pool = Pool(processes=6)
        for chr_tbp in obj.chrs_tbp:
            #compute_worker(obj, chr_tbp)
            #pool.apply(compute_worker, (obj, chr_tbp))
            pool.apply_async(compute_worker, (obj, chr_tbp))
        pool.close()
        pool.join()
         
        obj.file_combine()

def compute_worker(obj, chr_tbp):
    print chr_tbp
    h5 = tb.openFile(obj.h5)
    sample1 = h5.getNode("/", obj.track1)
#    sample_chrs = [chr._v_name for chr in sample1._f_iterNodes()]
    #if not chr_tbp in sample_chrs: return
    sample2 = h5.getNode("/", obj.track2)
    sample1_data = sample1._f_getChild(chr_tbp)
    sample2_data = sample2._f_getChild(chr_tbp)
    track_wsize = int(sample1_data.getAttr('Resolution'))
    #pdb.set_trace()
    out_wig = open(obj.tmp_path + "/" + chr_tbp, 'w')
    out_wig.write("fixedStep chrom={0} start=1 step={1} span={2}\n".format(chr_tbp, obj.step, obj.window))
    
    window = obj.window / track_wsize
    step = obj.step / track_wsize
    start = 0
    end = window
    vals1 = 0
    vals2 = 0
    corr = 0
    pvalue = 0
    
    sample1_data = sample1_data[:]
    sample2_data = sample2_data[:]
    
    while end < len(sample1_data):
        #pdb.set_trace()
        vals1 = sample1_data[start:end]
        vals2 = sample2_data[start:end]
        (corr, pvalue) = corr_wrapper(obj.func, vals1, vals2)
        #if math.isnan(corr): corr = 0
        out_wig.write(str(corr) + "\n")
        start = start + step
        end = start + window
        
    h5.close()
    out_wig.close

def corr_wrapper(func, *args):
    return(func(*args))
  
def val_ratio(v1, v2):

    
    val = float(sum(v1)) / sum(v2)
    #if val < 1: pdb.set_trace()
    #val = np.mean(ratio)
    if val == 0 or not np.isfinite(val): return((0,0))
    val = math.log(val, 2)
    #if val < 0 : pdb.set_trace()
    #if np.isnan(val) or not np.isfinite(val): val = 0
    #if not np.isfinite(val): pdb.set_trace()
    return((val,0))
    
def val_diff(v1, v2):
    return((np.mean(v1-v2),0))
    
def val_fraction(v1, v2):
    v1_s = sum(v1)
    v2_s = sum(v2)
    if v1_s < 10 and v2_s < 10: return(("NaN", 0))
    val = v1_s / (v1_s + sum(v2))
    #if np.isnan(val): val = "NA"
    return((val,0))
    
def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument(dest="file")
    parser.add_argument(dest="track1")
    parser.add_argument(dest="track2")
    parser.add_argument('-w', dest="window", default="1000")
    parser.add_argument('-s', dest="step", default="500")
    parser.add_argument("--data_type")
    parser.add_argument("--method", choices=['pearson', 'spearman', 'ratio', 'diff', 'fraction'])
    
    args = parser.parse_args()
        
    obj = corr_genome(args.file, args.track1, args.track2, args.data_type, args.method, args.window, args.step)
    compute_by_chrom(obj)

    

if __name__ == "__main__":
    main(sys.argv)
