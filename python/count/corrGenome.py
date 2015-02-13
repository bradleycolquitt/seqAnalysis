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

out_path = "/media/storage2/analysis/genomecor"

## Correlate two tracks in sliding windows across genome
## Input h5, two tracks, window size, step size
## Output as wig

class corr_genome:
    def __init__(self, h5, track, data_type, func, window, step):
        self.h5 = h5
        self.track = track

        self.func = ""
        if func == "energy": self.func = val_energy
        #if method == "spearman": self.func = ss.spearmanr
        self.window = int(window)
        self.step = int(step)
        self.tmp_path = tempfile.mkdtemp()
        self.out_path = "/".join([out_path, data_type, func])
        if not os.path.exists(self.out_path): os.makedirs(self.out_path)
        self.outfile = "/".join([self.out_path, "_".join([track, "W" + window, "S" + step]) + ".wig"])
        if os.path.exists(self.outfile):
            dec = raw_input("File exists. Overwrite? [y/n]")
            if dec == "n": sys.exit()
        
        # get chrs tbp
        h5_data = tb.openFile(h5)
        sample1 = h5_data.getNode("/", track)
#        sample2 = h5_data.getNode("/", track2)
        self.chrs_tbp = [chr._v_name for chr in sample1._f_iterNodes()]
 #       sample2_chrs = [chr._v_name for chr in sample2._f_iterNodes()]
        #pdb.set_trace()
  #      self.chrs_tbp = list(set(sample1_chrs) & set(sample2_chrs))
        
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
    sample1 = h5.getNode("/", obj.track)
#    sample_chrs = [chr._v_name for chr in sample1._f_iterNodes()]
    #if not chr_tbp in sample_chrs: return
    #sample2 = h5.getNode("/", obj.track2)
    sample1_data = sample1._f_getChild(chr_tbp)
    #sample2_data = sample2._f_getChild(chr_tbp)
    track_wsize = int(sample1_data.getAttr('Resolution'))
    #pdb.set_trace()
    out_wig = open(obj.tmp_path + "/" + chr_tbp, 'w')
    out_wig.write("fixedStep chrom={0} start=1 step={1} span={2}\n".format(chr_tbp, obj.step, obj.window))
    
    window = obj.window / track_wsize
    step = obj.step / track_wsize
    start = 0
    end = window
    vals1 = 0
    #vals2 = 0
    corr = 0
    pvalue = 0
    
    
    while end < len(sample1_data):   
        #pdb.set_trace()
        vals1 = sample1_data[start:end]
        #vals2 = sample2_data[start:end]
        #(corr, pvalue) = corr_wrapper(obj.func, vals1, vals2)
        #if math.isnan(corr): corr = 0
        out_wig.write(str(func_wrapper(obj.func, vals1)) + "\n")
        start = start + step
        end = start + window
        
    h5.close()
    out_wig.close

def func_wrapper(func, *args):
    return(func(*args))

def val_mean(vals):
    result = np.mean(vals)
    if np.isnan(result): result = 0
    return result

def val_energy(vals):
    m = val_mean(vals)
    result = sum(pow(vals - m, 2))
    return result
    
def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument(dest="file")
    parser.add_argument(dest="track")
    #parser.add_argument(dest="track2")
    parser.add_argument('-w', dest="window", default=1000)
    parser.add_argument('-s', dest="step", default=500)
    parser.add_argument("--data_type")
    parser.add_argument("--function")
    
    args = parser.parse_args()
        
    obj = corr_genome(args.file, args.track, args.data_type, args.function, args.window, args.step)
    compute_by_chrom(obj)

    

if __name__ == "__main__":
    main(sys.argv)
