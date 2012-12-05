#! /usr/bin/env python

import sys
import os
import shutil
import argparse
import tempfile
import pdb
import scipy.signal as ss
import scipy.stats.stats as stats
import tables as tb
from multiprocessing import Pool

feature_path = "/home/user/lib/features_general/"
out_path = "/media/storage2/analysis/crosscor"

class cross_track:
    
    def __init__(self, h5_a, h5_b, track1, track2, feature, data_type, corr_type, flank):
        self.h5_a = h5_a
        self.h5_b = h5_b
        self.track1 = track1
        self.track2 = track2
        self.feature = feature_path + feature
        self.flank = flank
        self.tmp_path = tempfile.mkdtemp(suffix=os.path.basename(feature))
        self.corr_type = corr_type
        if track2 == track1: self.corr_type = "auto"
        self.func = ss.fftconvolve
        if corr_type == "spearmanr": self.func = stats.spearmanr
        self.out_path = "/".join([out_path, data_type, corr_type, feature])
        if not os.path.exists(self.out_path): os.makedirs(self.out_path)
        self.outfile = "/".join([self.out_path, "_".join([track1, track2, str(flank)])])
        if os.path.exists(self.outfile):
            dec = raw_input("File exists. Overwrite? [y/n]")
            if dec == "n": sys.exit()
            
        self.chrs_tbp = os.listdir(self.feature)
    
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
        pool = Pool(processes=4)
        for chr_tbp in obj.chrs_tbp:
            compute_worker(obj, chr_tbp)
            #pool.apply(compute_worker, (obj, chr_tbp))
            #pool.apply_async(compute_worker, (obj, chr_tbp))
        pool.close()
        pool.join()
        
        obj.file_combine()

def compute_worker(obj, chr_tbp):
    print chr_tbp
    h5_a = tb.openFile(obj.h5_a)
    h5_b = tb.openFile(obj.h5_b)
    sample1 = h5_a.getNode("/", obj.track1)
    sample_chrs = [chr._v_name for chr in sample1._f_iterNodes()]
    if not chr_tbp in sample_chrs: return
    sample2 = h5_b.getNode("/", obj.track2)
    sample1_data = sample1._f_getChild(chr_tbp)
    sample2_data = sample2._f_getChild(chr_tbp)

    feature_data = open(obj.feature + "/" + chr_tbp)
    feature_out_path = obj.tmp_path + "/" + chr_tbp
    feature_out = open(feature_out_path, 'w')
   
    start = 0
    end = 0
    vals1 = 0
    vals2 = 0
    corr = 0
    for line in feature_data:
        line = line.strip()
        sline = line.split()
        #pdb.set_trace()
        start = int(sline[1]) - 1 - obj.flank
        end = int(sline[2]) + obj.flank
        vals1 = sample1_data[start:end]
        vals2 = sample2_data[start:end]
        
        if obj.corr_type == "cross" or obj.corr_type == "auto": 
            corr = ss.fftconvolve(vals1, vals2, 'same')
        elif obj.corr_type == "spearmanr":
            corr = [stats.spearmanr(vals1, vals2)[0]]
        elif obj.corr_type == "pearsonr":
            corr = [stats.pearsonr(vals1, vals2)[0]]
        feature_out.write("\t".join(sline[3:5] + map(str,corr)) + "\n")
       
    feature_data.close()
    feature_out.close()
    h5_a.close()
    h5_b.close()
    
def corr_wrapper(func, *args):
    return(func(*args))
    
def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument(dest="feature")
    parser.add_argument(dest="file")
    parser.add_argument(dest="track1")
    parser.add_argument("--file_b", required=False)
    parser.add_argument("-b", dest="track2", required=False)
    parser.add_argument("--data_type")
    parser.add_argument("--cor_type", choices=['cross', 'spearmanr', 'pearsonr'])
    parser.add_argument("--flank", dest="flank", type=int, required=False, default=0)
    
    args = parser.parse_args()
    
    file_b = args.file
    if args.file_b: file_b = args.file_b
    
    track2 = ""
    if args.track2:
        # cross-correlate
        track2 = args.track2
    else:
        # auto-correlate
        track2 = args.track1
        
    obj = cross_track(args.file, file_b, args.track1, track2, args.feature, args.data_type, args.cor_type, args.flank)
    compute_by_chrom(obj)

    

if __name__ == "__main__":
    main(sys.argv)
