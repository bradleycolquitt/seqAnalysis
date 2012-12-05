#! /usr/bin/env python

import sys
import os
import argparse
import re
import tempfile
import shutil
import numpy as np
from scipy import stats
from string import atoi, atof
from multiprocessing import Pool

FEATURE_PATH = '/home/user/lib/features_general'
SAMPLE_PATH = '/media/storage2/data/wig'
OUT_PATH = '/media/storage2/analysis/features'
SAMPLES = []
"""
CELLS_SAMPLES = ['omp_hmedip.bed', 'ngn_hmedip.bed', 'icam_hmedip.bed',
                 'omp_medip.bed', 'ngn_medip.bed', 'icam_medip.bed']
"""
CELLS_SAMPLES = ['omp_hmc_rlm', 'ngn_hmc_rlm', 'icam_hmc_rlm',
                 'omp_mc_rlm', 'ngn_mc_rlm', 'icam_mc_rlm']
"""
D3A_SAMPLES = ['moe_wt_mc.bed', 'moe_d3a_mc.bed',
               'moe_wt_hmc.bed', 'moe_d3a_hmc.bed']
"""
D3A_SAMPLES = ['moe_wt_mc_rlm', 'moe_d3a_mc_rlm',
               'moe_wt_hmc_rlm', 'moe_d3a_hmc_rlm']


class comparison:
    def __init__(self, feature, date, samples, window, fun):
        #self.feature_path = feature_paths
        #self.feature_files = [f for f in os.listdir(self.feature_path) if re.search("chr", f)]
        self.feature = "/".join([FEATURE_PATH, feature])
        self.date = date
        self.sample_paths = ["/".join([SAMPLE_PATH, date, sample, window]) for sample in samples]
        self.sample_path2 = "/".join([SAMPLE_PATH, date, samples[1], window])
        self.samples = samples
        self.sample2 = samples[1]
        self.window = window
        self.fun = fun
        self.out_path = "/".join([OUT_PATH, "_".join([self.sample1, self.sample2]), os.path.basename(feature)])
        if os.path.exists(self.out_path):
            print "file exists: " + self.out_path
            return
       
    def run(self):
        feature = self.feature
        samples = self.sample_paths
        feature_chrs = os.listdir(feature)
        samples_chrs = [os.listdir(sample) for sample in samples]
        chrs_tbp = [list(set(feature_chrs) & set(sample_chrs)) for sample_chrs in samples_chrs]
        tmp_path = tempfile.mkdtemp(suffix=os.path.basename(feature))
        
        for chr_tbp in chrs_tbp:
            feature_data = open(feature + "/" + chr_tbp)
            sample_data = open(sample + "/" + chr_tbp)
            feature_out_path = tmp_path + "/" + chr_tbp
            if os.path.exists(feature_out_path): 
                print "continue"
                continue
            feature_out = open(feature_out_path, 'w')
            self.tab_core(feature_data, sample_data, feature_out, fun)
            feature_data.close()
            sample_data.close()
            feature_out.close()    
        self.file_combine(tmp_path)

    def tab_core(self, feature_data, sample_data, feature_out):
        window_size = atoi(self.window)
        feature_line = feature_data.readline().strip().split()
            
        ## Read in WIG file
        ## Read in feature file
        ## Generate WIG list with:
        ## Grab WIG value by each range of feature
        val_list = []
        for line in sample_data:
            if re.search("Step", line): continue
            val_list.append(atof(line.strip()))
        #print type(val_list[0])
        for line in feature_data:
            #if re.search("Step", line): continue
            line = line.split()
            start = atoi(line[1]) / window_size
            end = atoi(line[2]) / window_size
            vals = val_list[start:end]
            result = 0
            if fun == "mean":
                result = np.mean(vals)
            
               
            out = "\t".join([line[0], line[1], line[2], line[3], str(result)]) + "\n"
            feature_out.write(out)
        feature_out.close()
        
    def file_combine(self, tmp_path):
        feature = self.feature
        sample = self.sample
        out_prefix = OUT_PATH + "/" + os.path.basename(sample)
        print out_prefix
        if not os.path.exists(out_prefix): os.makedirs(out_prefix)
        out = open(out_prefix + "/" + os.path.basename(feature), 'wa')
        #tmp_path = sample + "/" + os.path.basename(feature) + "_tmp"
        files_tmp = os.listdir(tmp_path)
        for file_tmp in files_tmp:
            a = open(tmp_path + "/" + file_tmp)
            for line in a:
                out.write(line)
        shutil.rmtree(tmp_path)
        out.close()

    
    
def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', dest='samples', nargs=2, required=False, default="")
    parser.add_argument('--feature_set', dest='feature_set', required=False)
    parser.add_argument('-f', dest='feature', required=False)
    parser.add_argument('--sample_set', dest='sample_set', required=False)
    parser.add_argument('-d', dest='date')
    parser.add_argument('-w', dest='window')
    parser.add_argument('--function', dest='fun', required=False)
    args = parser.parse_args()
    
    feature_path = '/home/user/lib/features_general/'
    feature = args.feature
    fun = ""
    if not args.fun:
        fun = "mean"
    else:
        fun = args.fun
    #if args.sample_set:
    #    feature_files = [f for f in os.listdir(FEATURE_PATH) if re.search("chr", f)]
    #    samples = [args.sample]
    #    if args.sample_set == "cells":
    #        samples = CELLS_SAMPLES
    #    elif args.sample_set == "d3a":
    #        samples = D3A_SAMPLES
    #    pool = Pool(processes=len(samples))
    #    for feature in feature_files:
    #        result = [pool.apply_async(tab_worker, (feature, args.date, sample, args.window, args.fun)) \
    #                  for sample in samples]
    #        pool.close()
    #        pool.join()
        
    obj = comparison(args.feature, args.date, args.samples, args.window, args.fun)
    #obj = tab(feature_path, feature, samples, args.window) 
    obj.run()

   
if __name__ == "__main__":
    main(sys.argv)
