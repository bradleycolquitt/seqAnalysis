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
import numpy as np
from scipy import stats
from string import atoi, atof
from multiprocessing import Pool

ANNO_PATH = '/home/user/lib/annotations'
FEATURE_PATH = '/home/user/lib/features_general'
SAMPLE_PATH = '/media/storage2/data/h5'
ANNO_OUT_PATH = '/media/storage2/analysis/profiles/norm'
FEATURE_OUT_PATH = '/media/storage2/analysis/features/norm'
SAMPLES = []
"""
CELLS_SAMPLES = ['omp_hmedip.bed', 'ngn_hmedip.bed', 'icam_hmedip.bed',
                 'omp_medip.bed', 'ngn_medip.bed', 'icam_medip.bed']
"""
#CELLS_SAMPLES = ['omp_hmc_rlm', 'ngn_hmc_rlm', 'icam_hmc_rlm',
#                 'omp_mc_rlm', 'ngn_mc_rlm', 'icam_mc_rlm']

"""
D3A_SAMPLES = ['moe_wt_mc.bed', 'moe_d3a_mc.bed',
               'moe_wt_hmc.bed', 'moe_d3a_hmc.bed']
"""
D3A_SAMPLES = ['moe_wt_mc_rlm', 'moe_d3a_mc_rlm',
               'moe_wt_hmc_rlm', 'moe_d3a_hmc_rlm']


class tab:
    def __init__(self, anno, h5, fun, type, split_anno):
        self.anno = anno
        self.h5 = h5
        self.fun = fun
        self._type = type
        self.split_anno = split_anno
        self.out_path = ""
        if type == "anno":
            self.out_path = "/".join([ANNO_OUT_PATH, os.path.basename(anno)])
        elif type == "feature":
            self.out_path = "/".join([FEATURE_OUT_PATH, os.path.basename(anno)])
        if split_anno:
            self.out_path = self.out_path + "_split"
        if not os.path.exists(self.out_path): os.makedirs(self.out_path)
        
    def run(self):
        anno = self.anno
        for sample in self.h5.iterNodes("/"):
            if os.path.exists(self.out_path + "/" + sample._v_name):
                print "File exists. Skipping..."
                continue
            print sample._v_name
            anno_chrs = os.listdir(anno)
            sample_chrs = [chr._v_name for chr in sample._f_iterNodes()]
            chrs_tbp = list(set(anno_chrs) & set(sample_chrs))
            tmp_path = tempfile.mkdtemp(suffix=os.path.basename(anno))
            #worker(anno, sample, chrs_tbp, tmp_path, self.fun)
            #pool = Pool(processes=6)
            for chr_tbp in chrs_tbp:
                print chr_tbp
                
                self.tab_core(anno, sample, chr_tbp, tmp_path, self.fun)
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

    def tab_core(self, anno, sample, chr_tbp, tmp_path, fun):
        anno_data = open(anno + "/" + chr_tbp)
        sample_data = sample._f_getChild(chr_tbp)
        anno_out_path = tmp_path + "/" + chr_tbp
        anno_out = open(anno_out_path, 'w')
        window_size = int(sample_data.getAttr('Resolution'))
        anno_line = anno_data.readline().strip().split()
        #pdb.set_trace()    
        for line in anno_data:
            line = line.strip()
            sline = line.split()
            start = atoi(sline[1]) / window_size
            end = atoi(sline[2]) / window_size
            vals = sample_data[start:end]
            if len(vals) > 0:
                if self.split_anno:
                    for val in vals:
                        out = "\t".join([line, str(val)]) + "\n"
                        anno_out.write(out) 
                else:
                    vals = vals[vals>0]
                    if len(vals) > 0:
                        result = np.mean(vals)
                    #if np.isnan(result): result = 0
                        out = "\t".join([line, str(result)]) + "\n"
                        anno_out.write(out)
        anno_data.close()
        anno_out.close()
            
    def file_combine(self, tmp_path, sample_name):
        anno = self.anno
        out = open("/".join([self.out_path, sample_name]), 'w')
        files_tmp = os.listdir(tmp_path)
        for file_tmp in files_tmp:
            a = open(tmp_path + "/" + file_tmp)
            for line in a:
                out.write(line)
        shutil.rmtree(tmp_path)
        out.close()
    
def worker(anno, sample, chrs_tbp, tmp_path, fun):
    pool = Pool(processes=6)
    for chr_tbp in chrs_tbp:
        print chr_tbp
        #tab_core(anno, sample, chr_tbp, tmp_path, self.fun)
        pool.apply_async(tab_core, (anno, sample, chr_tbp, tmp_path, fun))
    pool.close()
    pool.join()
                    
                #self.tab_core(anno_data, sample_data, self.fun)
                #anno_data.close()
                #sample_data.close()
                #anno_out.close()    
    file_combine(anno, tmp_path, sample._v_name)
def tab_worker(anno, h5, fun, type, split_anno):
    print anno
    obj = tab(anno, h5, fun, type, split_anno)
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
        
    h5 = tb.openFile(SAMPLE_PATH + "/"+ args.sample_h5, 'r')
    
    if args.anno_set:
        annos = [ANNO_PATH + "/"+ f for f in os.listdir(ANNO_PATH) if re.search("chr", f)]
        [tab_worker(anno, h5, fun, "anno", args.split_anno) for anno in annos]
        #pool = Pool(processes=2)
        #result = [pool.apply_async(tab_worker, (anno, h5, args.fun)) \
        #        for anno in annos]
        #pool.close()
        #pool.join()
    elif args.feature_set:
        print args.split_anno
        features = [FEATURE_PATH + "/" + f for f in os.listdir(FEATURE_PATH) if re.search("chr", f)]
        print features
        [tab_worker(feature, h5, fun, "feature", args.split_anno) for feature in features]
    elif args.annotation:       
        obj = tab(ANNO_PATH + "/" + args.annotation, h5, fun, "anno", args.split_anno)
        obj.run()
    elif args.feature:
        obj = tab(FEATURE_PATH + "/" + args.feature, h5, fun, "feature", args.split_anno)
        obj.run()

   
if __name__ == "__main__":
    main(sys.argv)
