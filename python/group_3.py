#! /usr/bin/env python

import sys
import os
import argparse
import re
import tables as tb
import tempfile
import shutil
import numpy as np
from scipy import stats
from string import atoi, atof
from multiprocessing import Pool

ANNO_PATH = '/home/user/lib/annotations'
FEATURE_PATH = '/home/user/lib/features'
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


def tab(anno, h5, fun, type):
    anno = "/".join([ANNO_PATH, anno])
    out_path = ""
    if type == "anno":
        out_path = "/".join([ANNO_OUT_PATH, os.path.basename(anno)])
    elif type == "feature":
        out_path = "/".join([FEATURE_OUT_PATH, os.path.basename(anno)])
    if not os.path.exists(out_path): os.makedirs(out_path)
    for sample in h5.iterNodes("/"):
        sample_name = sample._v_name
        if os.path.exists(out_path + "/" + sample_name): continue
        print sample_name
        anno_chrs = os.listdir(anno)
        sample_chrs = [chr._v_name for chr in sample._f_iterNodes()]
        chrs_tbp = list(set(anno_chrs) & set(sample_chrs))
        print anno_chrs, sample_chrs, chrs_tbp
        tmp_path = tempfile.mkdtemp(suffix=os.path.basename(anno))
        pool = Pool(processes = 6)
        [pool.apply_async(tab_core, (anno, sample, chr, tmp_path, fun)) for chr in chrs_tbp]
        pool.close()
        pool.join()
        file_combine(tmp_path, out_path, sample_name)

def tab_core(anno, sample, chr_tbp, tmp_path, fun):
    anno_data = open(anno + "/" + chr_tbp)
    sample_data = sample._f_getChild(chr_tbp)
    anno_out_path = tmp_path + "/" + chr_tbp
    anno_out = open(anno_out_path, 'w')
    window_size = int(sample_data.getAttr('Resolution'))
    anno_line = anno_data.readline().strip().split()
        
    for line in anno_data:
        line = line.strip()
        sline = line.split()
        start = atoi(sline[1]) / window_size
        end = atoi(sline[2]) / window_size
        vals = sample_data[start:end]
        result = np.mean(vals)
        out = "\t".join([line, str(result)]) + "\n"
        anno_out.write(out)
    anno_data.close()
    anno_out.close()
        
def file_combine(tmp_path, out_path, sample_name):
    out = open("/".join([out_path, sample_name]), 'w')
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
def tab_worker(anno, h5, fun, type):
    print anno
    obj = tab(anno, h5, fun, type)
    obj.run()
    
def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', dest='annotation', required=False, default="")
    parser.add_argument('--anno_set', action="store_true", required=False)
    parser.add_argument('--feature_set', action="store_true", required=False)
    parser.add_argument('-t', dest='sample_h5')
    parser.add_argument('--function', dest='fun', required=False)
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
        [tab(anno, h5, fun, "anno") for anno in annos]
        #pool = Pool(processes=2)
        #result = [pool.apply_async(tab_worker, (anno, h5, args.fun)) \
        #        for anno in annos]
        #pool.close()
        #pool.join()
    elif args.feature_set:
        features = [FEATURE_PATH + "/" + f for f in os.listdir(FEATURE_PATH) if re.search("chr", f)]
        [tab(feature, h5, fun, "feature") for feature in features]
    else:
        tab(args.annotation, h5, args.fun, "anno")
        

   
if __name__ == "__main__":
    main(sys.argv)
