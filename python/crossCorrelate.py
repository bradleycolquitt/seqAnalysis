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
import numpy as np
import scipy.signal as ss
import scipy.stats.stats as stats
import tables as tb
from multiprocessing import Pool

feature_path = "/home/user/lib/features_general/"
out_path = "/media/storage2/analysis/crosscor"

class cross_track:
    
    def __init__(self, file_a, file_b, track_a, track_b, feature, data_type, corr_type, flank):
        self.file_a = file_a
        self.file_b = file_b
        self.track_a = track_a
        self.track_b = track_b
        self.input_data_type = "h5"
        #pdb.set_trace()
        if re.search("bam", file_a):
            self.track_a = file_a.split(".bam")[0]
            self.track_b = file_b.split(".bam")[0]
            self.input_data_type = "bam"
        self.feature = feature_path + feature
        self.flank = flank
        self.tmp_path = tempfile.mkdtemp(suffix=os.path.basename(feature))
        self.corr_type = corr_type
        if self.track_b == self.track_a: self.corr_type = "auto"
        self.func = ss.fftconvolve
        if corr_type == "spearmanr": self.func = stats.spearmanr
        #pdb.set_trace()
        self.out_path = "/".join([out_path, data_type, corr_type, feature])
        if not os.path.exists(self.out_path): os.makedirs(self.out_path)
        self.outfile = "/".join([self.out_path,
                                 "_".join([os.path.basename(self.track_a),
                                           os.path.basename(self.track_b),
                                           str(flank)])])
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
            if obj.input_data_type == "h5":
                compute_worker_h5(obj, chr_tbp)
            elif obj.input_data_type == "bam":
                compute_worker_bam(obj, chr_tbp)
            #pool.apply(compute_worker, (obj, chr_tbp))
            #pool.apply_async(compute_worker, (obj, chr_tbp))
        pool.close()
        pool.join()
        
        obj.file_combine()

def compute_worker_h5(obj, chr_tbp):
    print chr_tbp
    h5_a = tb.openFile(obj.file_a)
    h5_b = tb.openFile(obj.file_b)
    sample1 = h5_a.getNode("/", obj.track_a)
    sample_chrs = [chr._v_name for chr in sample1._f_iterNodes()]
    if not chr_tbp in sample_chrs: return
    sample2 = h5_b.getNode("/", obj.track_b)
    sample1_data = sample1._f_getChild(chr_tbp)
    sample2_data = sample2._f_getChild(chr_tbp)

    feature_data = open(obj.feature + "/" + chr_tbp)
    feature_out_path = obj.tmp_path + "/" + chr_tbp
    feature_out = open(feature_out_path, 'w')
   
    window_size = int(sample1_data.getAttr('Resolution'))
   
    start = 0
    end = 0
    vals1 = 0
    vals2 = 0
    corr = 0
    for line in feature_data:
        line = line.strip()
        sline = line.split()
        if re.search("chr3-438", line): pdb.set_trace()
        start = (int(sline[1]) - 1 - obj.flank) / window_size
        end = (int(sline[2]) + obj.flank) / window_size
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

def compute_worker_bam(obj, chr_tbp):
    print chr_tbp
    
    file_a = pysam.Samfile(obj.file_a, 'rb')
    file_b = pysam.Samfile(obj.file_b, 'rb')
    
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
        
        vals1 = np.zeros(end - start)
        vals2 = np.zeros(end - start)
      
        for column in file_a.pileup(reference=chr_tbp, start=start, end=end):
            if (column.pos >= start and column.pos < end):
 #              pdb.set_trace()
                try:
                    vals1[(column.pos - start)] = column.n
                except:
                    pdb.set_trace()
        
        for column in file_b.pileup(reference=chr_tbp, start=start, end=end):
            if (column.pos >= start and column.pos < end):
#               pdb.set_trace()
                try:
                    vals2[(column.pos - start)] = column.n
                except:
                    pdb.set_trace()            
        
        
        
        if obj.corr_type == "cross" or obj.corr_type == "auto": 
            corr = ss.fftconvolve(vals1, vals2, 'same')
        elif obj.corr_type == "spearmanr":
            corr = [stats.spearmanr(vals1, vals2)[0]]
        elif obj.corr_type == "pearsonr":
            corr = [stats.pearsonr(vals1, vals2)[0]]
        feature_out.write("\t".join(sline[3:5] + map(str,corr)) + "\n")
       
    feature_data.close()
    feature_out.close()
   

    
def corr_wrapper(func, *args):
    return(func(*args))
    
def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument(dest="feature", help=" ".join(["Feature in ", feature_path]))
    parser.add_argument("--file_a", dest="file_a", required=False)
    parser.add_argument("--track_a", dest="track_a", required=False)
    parser.add_argument("--file_b", required=False)
    parser.add_argument("--track_b", dest="track_b", required=False)
    parser.add_argument("--data_type")
    parser.add_argument("--cor_type", choices=['cross', 'spearmanr', 'pearsonr'])
    parser.add_argument("--flank", dest="flank", type=int, required=False,
                        default=0, help="Extend each bed record by this amount")
    
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
        
    obj = cross_track(args.file_a, file_b, args.track_a, track_b, args.feature, args.data_type, args.cor_type, args.flank)
    compute_by_chrom(obj)

    

if __name__ == "__main__":
    main(sys.argv)
