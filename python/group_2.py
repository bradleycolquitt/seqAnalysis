#! /usr/bin/env python

import sys
import os
import argparse
import pdb
import re
import tables as tb
import tempfile
import shutil
import warnings
import operator
import numpy as np
from math import sqrt
from scipy import stats
from string import atoi, atof
from multiprocessing import Pool

#ANNO_PATH = '/home/user/lib/annotations_hires'
ANNO_PATH = '/home/user/lib/annotations'
FEATURE_PATH = '/home/user/lib/features_general'
SAMPLE_PATH = '/media/storage2/data/h5'
ANNO_OUT_PATH = '/media/storage2/analysis/profiles/norm'
FEATURE_OUT_PATH = '/media/storage2/analysis/features/norm'

class tab:
    def __init__(self, anno, h5, fun, type, data_type, split_anno):
        self.anno = anno
        self.h5 = h5
        self.fun = fun
        self._type = type
        self.split_anno = split_anno
        self.window_size = 0
        self.out_path = ""
        if type == "anno":
            self.out_path = "/".join([ANNO_OUT_PATH, data_type, fun, os.path.basename(anno)])
        elif type == "feature":
            if not split_anno:
                self.out_path = "/".join([FEATURE_OUT_PATH, data_type, fun, os.path.basename(anno)])
            else:
                self.out_path = "/".join([FEATURE_OUT_PATH, data_type, "split", os.path.basename(anno)])    
        if not os.path.exists(self.out_path): os.makedirs(self.out_path)
        
    def run(self):
        anno = self.anno
        #pdb.set_trace()
        for sample in self.h5.iterNodes("/"):
            if os.path.exists(self.out_path + "/" + sample._v_name):
                print "File exists. Skipping..."
                continue
            print sample._v_name
            anno_chrs = os.listdir(anno)
            sample_chrs = [chr._v_name for chr in sample._f_iterNodes()]
            chrs_tbp = list(set(anno_chrs) & set(sample_chrs))
            tmp_path = tempfile.mkdtemp(suffix=os.path.basename(anno))
            
            exp = 0
            if self.fun=="randc":
                try:
                    exp = sample._f_getAttr('ExpectedRPM')
                except AttributeError:
                    exp = computeExpected(sample)
                    sample._f_setAttr('ExpectedRPM', exp)
                    self.h5.flush()
                    print "Expected RPM:", exp
            #worker(anno, sample, chrs_tbp, tmp_path, self.fun)
            #pool = Pool(processes=6)
            for chr_tbp in chrs_tbp:
                print chr_tbp
                self.tab_core(anno, sample, chr_tbp, tmp_path, self.fun, exp)
            ##    pool.apply_async(tab_core, (anno, sample, chr_tbp, self.fun))
            ##pool.close()
            ##pool.join()
            #    
            #    
            #    
            #    #self.tab_core(anno_data, sample_data, self.fun)
            #    #anno_data.close()
            #    #sample_data.close()
            #    #anno_out.close()    
            self.file_combine(tmp_path, sample._v_name)
        #self.h5.flush()
        
    def tab_core(self, anno, sample, chr_tbp, tmp_path, fun, exp):
        anno_data = open(anno + "/" + chr_tbp)
        sample_data = sample._f_getChild(chr_tbp)
        anno_out_path = tmp_path + "/" + chr_tbp
        anno_out = open(anno_out_path, 'w')
        self.window_size = int(sample_data.getAttr('Resolution'))
        anno_line = anno_data.readline().strip().split()
        #print fun    
        for line in anno_data:
            line = line.strip()
            sline = line.split()
            start = atoi(sline[1]) / self.window_size
            end = atoi(sline[2]) / self.window_size
            #pdb.set_trace()
            vals = sample_data[start:(end+1)]
            if len(vals) > 0:
                #pdb.set_trace()
                if self.split_anno:
                    for val in vals:
                        out = "\t".join([line, str(val)]) + "\n"
                        anno_out.write(out) 
                else:
                    result = 0
                    if fun == "mean":
                        result = val_mean(vals)
                    elif fun == "median":
                        result = np.median(vals)
                    elif fun == "mode":
                        result = val_mode(vals)
                    elif fun == "sum":
                        result = sum(vals)
                        #print result
                    elif fun == "var":
                        result = val_var(vals)
                    elif fun == "cv":
                        result = val_cv(vals)
                    elif fun == "randc":
                        result = val_compareRandom(vals, exp)
                    #if result > 4: pdb.set_trace()
                    out = "\t".join([line, str(result)]) + "\n"
                    anno_out.write(out)
                
        anno_data.close()
        anno_out.close()
            
    def file_combine(self, tmp_path, sample_name):
        anno = self.anno
        if self.split_anno:
            out = open("/".join([self.out_path, sample_name + "_" + str(self.window_size)]), 'w')
        else:
            out = open("/".join([self.out_path, sample_name]), 'w')
        files_tmp = os.listdir(tmp_path)
        for file_tmp in files_tmp:
            a = open(tmp_path + "/" + file_tmp)
            for line in a:
                out.write(line)
        shutil.rmtree(tmp_path)
        out.close()
        
def computeExpected(sample):
    print 'Computing expected RPM per window...'
    read_sum = 0
    window_num = 0
    for chrom in sample._f_iterNodes():
        read_sum = read_sum + sum(chrom)
        window_num = window_num + len(chrom)
    return(read_sum / window_num)
    
def val_mean(vals):
    result = np.mean(vals)
    if np.isnan(result): result = 0
    return result

def val_mode(vals):
    counts = {}
    for val in vals:
        if val in counts:
            counts[val] = counts[val] + 1
        else:
            counts[val] = 1
    return(max(counts.iteritems(), key=operator.itemgetter(1))[0])

def val_var(vals):
    mean = val_mean(vals)
    denom = 1
    if len(vals) > 1: denom = len(vals) - 1
    var = sum(pow(vals - mean, 2)) / denom 
    return(var)
    
def val_cv(vals):
    mean = val_mean(vals)
    if mean == 0: return(0)
    sd = sqrt(val_var(vals))
    return(sd / mean)
    
def val_compareRandom(vals, exp):
    val_sum = np.sum(vals)
    exp_sum = len(vals) * exp
    return(val_sum / exp_sum)
    
def worker(anno, sample, chrs_tbp, tmp_path, fun):
    pool = Pool(processes=6)
    for chr_tbp in chrs_tbp:
        print chr_tbp
        #tab_core(anno, sample, chr_tbp, tmp_path, self.fun)
        pool.apply_async(tab_core, (anno, sample, chr_tbp, tmp_path, fun))
    pool.close()
    pool.join()
    file_combine(anno, tmp_path, sample._v_name)
    
def tab_worker(anno, h5, fun, type, data_type, split_anno):
    print anno
    #db.set_trace()
    obj = tab(anno, h5, fun, type, data_type, split_anno)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj.run()
    
def main(argv):
    #warnings.simplefilter("ignore")
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', dest='annotation', required=False, default="")
    parser.add_argument('-f', dest='feature', required=False)
    parser.add_argument('--anno_set', action="store_true", required=False)
    parser.add_argument('--feature_set', action="store_true", required=False)
    parser.add_argument('-t', dest='sample_h5', required=False)
    parser.add_argument('--function', dest='fun', required=False)
    parser.add_argument('--split', dest="split_anno", action="store_true", default=False)
    parser.add_argument('--data_type')
    args = parser.parse_args()
    
    #annos = args.annotation
    fun = ""
    if not args.fun:
        fun = "mean"
    else:
        fun = args.fun
    
    h5 = tb.openFile(SAMPLE_PATH + "/"+ args.sample_h5, 'a')
    
    if args.anno_set:
        annos = [ANNO_PATH + "/"+ f for f in os.listdir(ANNO_PATH) if re.search("chr", f)]
        [tab_worker(anno, h5, fun, "anno", args.data_type, args.split_anno) for anno in annos]
        #pool = Pool(processes=2)
        #result = [pool.apply_async(tab_worker, (anno, h5, args.fun)) \
        #        for anno in annos]
        #pool.close()
        #pool.join()
    elif args.feature_set:
        
        features = [FEATURE_PATH + "/" + f for f in os.listdir(FEATURE_PATH) if re.search("chr", f)]
        #print features
        [tab_worker(feature, h5, fun, "feature", args.data_type, args.split_anno) for feature in features]
    elif args.annotation:
        #pdb.set_trace()
        tab_worker(ANNO_PATH + "/" + args.annotation, h5, fun, "anno", args.data_type, args.split_anno)
        #obj = tab(ANNO_PATH + "/" + args.annotation, h5, fun, "anno", args.data_type, args.split_anno)
        
        #obj.run()
    elif args.feature:
        tab_worker(FEATURE_PATH + "/" + args.feature, h5, fun, "feature", args.data_type, args.split_anno)
        #obj = tab(FEATURE_PATH + "/" + args.feature, h5, fun, "feature", args.data_type, args.split_anno)
        #obj.run()
   
if __name__ == "__main__":
    main(sys.argv)
