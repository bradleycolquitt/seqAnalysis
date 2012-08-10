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
import pysam
import numpy as np
from math import sqrt
from scipy import stats
from string import atoi, atof
from multiprocessing import Pool

ANNO_PATH = '/home/user/lib/annotations_hires'
#ANNO_PATH = '/home/user/lib/annotations'
FEATURE_PATH = '/home/user/lib/features_general'
SAMPLE_PATH = '/media/storage2/data/h5'
ANNO_OUT_PATH = '/media/storage2/analysis/profiles/norm'
FEATURE_OUT_PATH = '/media/storage2/analysis/features/norm'

class tab:
    def __init__(self, anno, sample, sample_type, fun, type, flank, data_type, split_anno):
        self.anno = anno
        self.sample = sample
        self.sample_type = sample_type
        self.fun = fun
        #self._type = type
        self.split_anno = split_anno
        self.window_size = 1
        self.flank = flank
        self.out_path = ""
        if type == "anno":
            self.out_path = "/".join([ANNO_OUT_PATH, data_type, fun, os.path.basename(anno)])
        elif type == "feature":
            if not split_anno:
                self.out_path = "/".join([FEATURE_OUT_PATH, data_type, fun, os.path.basename(anno)])
            else:
                self.out_path = "/".join([FEATURE_OUT_PATH, data_type, "split", os.path.basename(anno)])    
        if not os.path.exists(self.out_path): os.makedirs(self.out_path)
        
    def run_h5(self):
        anno = self.anno
        #pdb.set_trace()
        h5 = tb.openFile(SAMPLE_PATH + "/" + self.sample)
        samples = [s._v_name for s in h5.iterNodes("/")]
        h5.close()
        #for sample in self.sample.iterNodes("/"):
        for sample in samples:
            if self.flank:
                #if os.path.exists(self.out_path + "/" + sample._v_name + "_flank" + str(self.flank)):
                if os.path.exists(self.out_path + "/" + sample + "_flank" + str(self.flank)):
                    print "File exists. Skipping..."
                    continue
            else:
                #if os.path.exists(self.out_path + "/" + sample._v_name):
                if os.path.exists(self.out_path + "/" + sample):
                   print "File exists. Skipping..."
                   continue
            #print sample._v_name
            print sample
            #anno_chrs = os.listdir(anno)
            #sample_chrs = [chr._v_name for chr in sample._f_iterNodes()]
            chrs_tbp = os.listdir(anno)
            #chrs_tbp = list(set(anno_chrs) & set(sample_chrs))
            tmp_path = tempfile.mkdtemp(suffix=os.path.basename(anno))
            
            
            exp = 0
            #if self.fun=="randc":
            #    try:
            #        exp = sample._f_getAttr('ExpectedRPM')
            #    except AttributeError:
            #        exp = computeExpected(sample)
            #        sample._f_setAttr('ExpectedRPM', exp)
            #        self.sample.flush()
            #        print "Expected RPM:", exp
            #worker(anno, sample, chrs_tbp, tmp_path, self.fun)
            pool = Pool(processes=6)
            for chr_tbp in chrs_tbp:
                #print chr_tbp
                #self.tab_h5(anno, sample, chr_tbp, tmp_path, self.fun, exp)
                #tab_h5(self, anno, sample, chr_tbp, tmp_path, self.fun, exp)
                #pdb.set_trace()
                #pool.apply(tab_h5, (self, anno, sample, chr_tbp, tmp_path, self.fun, exp))
                pool.apply_async(tab_h5, (self, anno, sample, chr_tbp, tmp_path, self.fun, exp))
            pool.close()
            pool.join()
            #    
            #    
            #    
            #    #self.tab_core(anno_data, sample_data, self.fun)
            #    #anno_data.close()
            #    #sample_data.close()
            #    #anno_out.close()
            if self.flank == 0:
                #self.file_combine(tmp_path, sample._v_name)
                self.file_combine(tmp_path, sample)
            else:
                #self.file_combine(tmp_path, sample._v_name + "_flank" + str(self.flank))
                self.file_combine(tmp_path, sample + "_flank" + str(self.flank))
            
        #self.sample.flush()
        
    def run_bam(self):
        anno = self.anno
        #pdb.set_trace()
        if os.path.exists(self.out_path + "/" + self.sample.split(".bam")[0]):
            print "File exists. Skipping..."
            return
        #print self.sample
        
        chrs_tbp = os.listdir(anno)
        tmp_path = tempfile.mkdtemp(suffix=os.path.basename(anno))
      
        pool = Pool(processes=5)
        #bam = pysam.Samfile(self.sample, 'rb')
        for chr_tbp in chrs_tbp:
            #print chr_tbp
            #self.tab_bam(anno, self.sample, chr_tbp, tmp_path, self.fun, exp)
            #pool.apply_async(tab_bam, (anno, self.sample, chr_tbp, tmp_path, self.fun, exp))
            pool.apply_async(tab_bam, (self, chr_tbp, tmp_path))
            #tab_bam(self, chr_tbp, tmp_path)
        pool.close()
        pool.join()
        #pdb.set_trace()
        if self.flank == 0:
            self.file_combine(tmp_path, os.path.basename(self.sample).split(".bam")[0])
        else:
            self.file_combine(tmp_path, os.path.basename(self.sample).split(".bam")[0] + "_flank" + str(self.flank))
    def file_combine(self, tmp_path, sample_name):
        #pdb.set_trace()
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

        
def tab_h5(self, anno, track, chr_tbp, tmp_path, fun, exp):
    #pdb.set_trace()
    print chr_tbp
    h5 = tb.openFile(SAMPLE_PATH + "/" + self.sample)
    sample = h5.getNode("/", track)
    anno_data = open(anno + "/" + chr_tbp)
    sample_data = sample._f_getChild(chr_tbp)
    anno_out_path = tmp_path + "/" + chr_tbp
    anno_out = open(anno_out_path, 'w')
    self.window_size = int(sample_data.getAttr('Resolution'))
    anno_line = anno_data.readline().strip().split()
    flank = self.flank
    start = 0
    end = 0
    vals = ""
    for line in anno_data:
        line = line.strip()
        sline = line.split()
        if flank == 0:
            start = (atoi(sline[1]) - 1)
            end = (atoi(sline[2]) - 1)
            vals = sample_data[start:(end+1)]
            #pdb.set_trace()
            
        else:
            #pdb.set_trace()
            start = [int(sline[1]) - 1 - flank, int(sline[2])]
            end = [int(sline[1]) - 2, int(sline[2]) + flank - 1]
            #iter = range(start[0], end[0]) + range(start[1], end[1])
            #vals = [sample_data[i] for i in iter]
            
            vals = np.append(sample_data[start[0]:(end[0] + 1)], sample_data[start[1]:(end[1] + 1)])
            #pdb.set_trace()
        
        #pdb.set_trace()
        #vals = sample_data[start:(end+1)]
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
    h5.close()
        
#def tab_bam(self, anno, sample, chr_tbp, tmp_path, fun, exp):
def tab_bam(obj, chr_tbp , tmp_path):
    #pdb.set_trace()
    print chr_tbp
   
    anno_data = open(obj.anno + "/" + chr_tbp)
    bam = pysam.Samfile(obj.sample, 'rb')
    norm_val = float(bam.mapped) / 1E6
    
    anno_out_path = tmp_path + "/" + chr_tbp
    anno_out = open(anno_out_path, 'w')
    
    flank = obj.flank
    window_size = obj.window_size
    anno_line = anno_data.readline().strip().split()
    
    vals = 0
    
    #pdb.set_trace()
    start = 0
    end = 0
    for line in anno_data:
        line = line.strip()
        sline = line.split()
        if obj.flank == 0:
            start = [(atoi(sline[1]) - 1)]
            end = [(atoi(sline[2]) - 1)]
            vals = np.zeros(end[0] - start[0] + 1)
        else:
            start = [int(sline[1]) - 1 - flank, int(sline[2])]
            end = [int(sline[1]) - 2, int(sline[2]) + flank - 1]
            vals = np.zeros(flank * 2)
        #vals = np.zeros(end - start + 1)
        #pdb.set_trace()
        j = 0
        for i in xrange(len(start)):
            it = bam.pileup(chr_tbp, start[i], end[i] + 1)
            try:
                for proxy in it:
                    pos = proxy.pos - start[i]
                    if pos < 0: continue
                    if proxy.pos == end[i] + 1: break
                    vals[j] = proxy.n
                    j = j + 1
            except IndexError:
                pdb.set_trace()
                
        if obj.split_anno:
            for val in vals:
                out = "\t".join([line, str(val / norm_val)]) + "\n"
                anno_out.write(out) 
        else:
            result = 0
            fun = obj.fun
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
            out = "\t".join([line, str(result / norm_val)]) + "\n"
            
            anno_out.write(out)
        
    #pdb.set_trace()
    anno_data.close()
    anno_out.close()
        
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
    
def tab_worker(anno, sample, sample_type, fun, type, flank, data_type, split_anno):
    #print anno
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if sample_type == "h5":
            obj = tab(anno, sample, sample_type, fun, type, flank, data_type, split_anno)
            obj.run_h5()
        elif sample_type == "bam":
            for bam in sample:
                obj = tab(anno, bam, sample_type, fun, type, flank, data_type, split_anno)
                obj.run_bam()
    
def main(argv):
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', dest='annotation', required=False, default="", help="Annotation file")
    parser.add_argument('-f', dest='feature', required=False, help="Feature file")
    parser.add_argument('--anno_set', action="store_true", required=False, help="Process all annotations")
    parser.add_argument('--feature_set', action="store_true", required=False, help="Process all features")
    parser.add_argument('-t', dest='h5', required=False, help="HDF5 track file")
    parser.add_argument('-b', dest='bam', required=False, nargs='*', help="BAM file")
    parser.add_argument('--data_type', help="Directory for analysis output")
    parser.add_argument('--function', dest='fun', required=False, default="mean",
                        choices=['mean', 'median', 'mode','sum', 'var', 'cv'],
                        help="Summary function to apply to read counts")
    parser.add_argument('--split', dest="split_anno", action="store_true", default=False)
    parser.add_argument('--flank', type=int, required=False, default=0, help="Size of regions flanking bed to compute values for") 
    args = parser.parse_args()
    
    fun = args.fun
    
    sample = 0
    sample_type = ""
    if args.h5:
        #sample = tb.openFile(SAMPLE_PATH + "/"+ args.h5, 'a')
        sample = args.h5
        sample_type = "h5"
    elif args.bam:
        sample = args.bam
        sample_type = "bam"
        
    if args.anno_set:
        annos = [ANNO_PATH + "/"+ f for f in os.listdir(ANNO_PATH) if re.search("chr", f)]
        [tab_worker(anno, sample, sample_type, fun, "anno", args.data_type, args.split_anno) for anno in annos]
        
    elif args.feature_set:
        features = [FEATURE_PATH + "/" + f for f in os.listdir(FEATURE_PATH) if re.search("chr", f)]
        [tab_worker(feature, sample, sample_type, fun, "feature", args.data_type, args.split_anno) for feature in features]
    
    elif args.annotation:
        tab_worker(ANNO_PATH + "/" + args.annotation, sample, sample_type, fun, "anno", args.flank, args.data_type, args.split_anno)
    
    elif args.feature:
        tab_worker(FEATURE_PATH + "/" + args.feature, sample, sample_type, fun, "feature", args.flank, args.data_type, args.split_anno)
    
    
   
if __name__ == "__main__":
    main(sys.argv)
